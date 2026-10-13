"""Model construction records using resolved anatomy-builder settings."""

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

from ogo.fea.spine import (
    DEFAULT_SPINE_ISO_RESOLUTION_MM,
    DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    DEFAULT_SPINE_PMMA_INTRUSION_MM,
    DEFAULT_SPINE_PMMA_THICKNESS_MM,
    DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
    DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
    DEFAULT_SPINE_REGISTRATION_BACKEND,
    DEFAULT_SPINE_REGISTRATION_ITERATIONS,
    DEFAULT_SPINE_REGISTRATION_LANDMARKS,
    DEFAULT_SPINE_ICP_TARGET,
    SPINE_ALIGNMENT_METHOD,
    default_spine_reference_path,
    generation_settings as spine_generation_settings,
    solve_report_profile as spine_solve_report_profile,
)
from ogo.fea.femur import (
    DEFAULT_FEMUR_CUT_MODE,
    DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM,
    DEFAULT_FEMUR_GREATER_TROCHANTER_INCLUSION_LENGTH_MM,
    DEFAULT_FEMUR_ISO_RESOLUTION_MM,
    DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    DEFAULT_FEMUR_REGISTRATION_BACKEND,
    DEFAULT_FEMUR_REGISTRATION_ITERATIONS,
    DEFAULT_FEMUR_REGISTRATION_LANDMARKS,
    DEFAULT_PMMA_INTRUSION_MM,
    DEFAULT_PMMA_THICKNESS_MM,
    POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION,
    FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
    GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION,
    SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION,
    solve_report_profile as femur_solve_report_profile,
)

def option_value(options, name, default=None):
    """Read an already-parsed builder setting, with material fallback support."""
    value = options.get(name.lstrip("-"))
    return default if value is None else value


def option_float(options, name, default):
    return float(option_value(options, name, default))


def option_int(options, name, default):
    return int(option_value(options, name, default))


def option_optional_float(options, name):
    value = option_value(options, name)
    return None if value is None else float(value)


def option_path(options, name):
    return option_value(options, name)


def percent_displacement_metadata(
    model_path: Path,
    *,
    report_profile: str,
    failure_axis: str,
    target_percent: float,
) -> dict:
    """Describe the solve displacement derived from model geometry."""
    metadata = {
        "value_mm": None,
        "target_displacement_percent": target_percent,
        "characteristic_length_mm": None,
        "value_source": "target_displacement_percent * characteristic_length_mm / 100",
    }
    try:
        from ogo.util.faim import (
            infer_profile_characteristic_length_mm,
            read_prescribed_displacement,
        )

        characteristic_length = infer_profile_characteristic_length_mm(
            model_path,
            report_profile,
            failure_axis,
        )
        current_value = read_prescribed_displacement(model_path, report_profile)
        sign = -1.0 if current_value != "" and float(current_value) < 0 else 1.0
        target_mm = abs(float(target_percent)) * characteristic_length / 100.0
        metadata["value_mm"] = sign * target_mm
        metadata["characteristic_length_mm"] = characteristic_length
    except Exception as exc:
        metadata["value_note"] = "available after model generation with netCDF4: {}".format(exc)
    return metadata


def absolute_displacement_metadata(options: dict, default_value: str) -> dict:
    """Describe a solve displacement kept as an explicit absolute mm value."""
    value = option_value(options, "--fe_displacement", default_value)
    try:
        value_mm = float(value)
    except (TypeError, ValueError):
        value_mm = value
    return {
        "value_mm": value_mm,
        "target_displacement_percent": None,
        "characteristic_length_mm": None,
        "value_source": "explicit --fe_displacement used as absolute model displacement",
    }



def write_modeling_metadata(
    model_path: Path,
    model_type: str,
    generator_argv: Sequence[str],
    args: argparse.Namespace,
    bc_audit_summary: Optional[dict] = None,
    *, reporting: dict,
) -> Optional[Path]:
    """Write a traceable model-building record next to a generated n88model."""
    if args.dry_run:
        return None
    if not model_path.exists():
        return None

    output_path = model_path.with_name(model_path.with_suffix("").name + "_modeling.json")
    generation_settings = {}
    if model_type == "spine":
        generation_settings = spine_generation_settings()
        if output_path.exists():
            generation_settings = json.loads(output_path.read_text()).get(
                "generation_settings", generation_settings
            )
    mask_source = args.pistoia_mask or (args.bone_mask if args.pistoia_mask_label else None)
    mask_output = model_path.with_name(model_path.stem + "_pistoia_mask.nii.gz")
    from ogo.fea import femur, spine
    builder = spine if model_type == "spine" else femur
    options = vars(builder.build_parser().parse_args(list(generator_argv)))
    common = {
        "schema_version": 1,
        "model_file": str(model_path),
        "generator": {
            "entry_point": "ogoFEA",
            "model_type": model_type,
            "lower_level_argv": list(generator_argv),
        },
        "inputs": {
            "calibrated_image": str(args.calibrated_image),
            "bone_mask": str(args.bone_mask),
            "output_path": str(args.output_path) if args.output_path is not None else str(args.calibrated_image.parent),
        },
        "geometry": {
            "model_coordinates": "preprocessed_image_physical_space",
            "origin_policy": "cropping and resampling update the image origin; FE meshing preserves that origin",
            "boundary_condition_coordinates": "defined from the generated model bounding box in the same physical space",
        },
        "post_generation_validation": {
            "bc_audit": {
                "enabled": bc_audit_summary is not None,
                "flat_tolerance": args.bc_audit_flat_tolerance,
                "summary": bc_audit_summary,
                "json_path": None,
                "csv_path": None,
                "png_path": None
                if not args.debug
                else str(model_path.with_name(model_path.with_suffix("").name + "_bc_audit.png")),
            }
        },
        "solve_and_reporting": {
            "solve_requested": not args.no_solve,
            "critical_volume_pct": reporting["critical_volume_pct"],
            "masked_critical_volume_pct": (
                reporting["critical_volume_pct"]
                if getattr(args, "masked_critical_volume", None) is None
                else args.masked_critical_volume
            ),
            "critical_strain": reporting["critical_strain"],
            "target_displacement_percent": reporting["target_displacement_percent"],
            "target_displacement_definition": "percent strain converted from model geometry for spine; percent of femur length for femur",
            "run_pistoia": args.run_pistoia or args.require_pistoia or mask_source is not None,
            "source_pistoia_mask": None
            if mask_source is None
            else str(mask_source),
            "pistoia_mask_labels": list(args.pistoia_mask_label or []),
            "model_space_pistoia_mask": None
            if mask_source is None
            else str(mask_output),
            "compress_solved_model": not args.no_compress,
        },
    }


    build = _spine_metadata if model_type == "spine" else _femur_metadata
    metadata = build(common, model_path, options, args, generation_settings)
    metadata["generation_arguments"] = options
    if generation_settings:
        metadata["generation_settings"] = generation_settings
    output_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {output_path}")
    return output_path


def material_metadata(options, site, *, include_cortical):
    """Record regional material laws and the shared PMMA settings."""
    trabecular = {
        "elastic_E_func": option_value(options, "--elastic_E_func", "default_E"),
        "yield_comp_func": option_value(options, "--yield_comp_func"),
        "yield_tens_func": option_value(options, "--yield_tens_func"),
        "material_id_range": [1, 128],
    }
    cortical = None
    if include_cortical:
        cortical = {
            key: option_value(options, "--cort_" + key, trabecular[key])
            for key in ("elastic_E_func", "yield_comp_func", "yield_tens_func")
        }
        cortical["poissons_ratio"] = option_float(
            options, "--cort_poissons_ratio", option_float(options, "--poissons_ratio", 0.3)
        )
        cortical["material_id_range"] = [129, 256]
    record = {
        "builder": "ogo.fea.materials.build_" + site + "_material_table",
        "convention": (
            "region material IDs; 0 background, 1..128 trabecular, "
            + ("optional " if site == "femur" else "")
            + "129..256 cortical, PMMA explicit"
        ),
        "poissons_ratio": option_float(options, "--poissons_ratio", 0.3),
        "trabecular": trabecular,
        "cortical": cortical,
        "pmma": {
            "material_id": option_int(options, "--pmma_mat_id", 5000),
            "elastic_E_MPa": option_float(options, "--pmma_E", 2500),
            "poissons_ratio": option_float(options, "--pmma_v", 0.3),
            "yield_compression_MPa": option_optional_float(options, "--pmma_yield_compression"),
            "yield_tension_MPa": option_optional_float(options, "--pmma_yield_tension"),
        },
    }
    if site == "femur":
        record["include_cortical_region"] = include_cortical
    return record


def _spine_metadata(common, model_path, options, args, generation_settings):
    body_label = option_int(options, "--mask_threshold", 0)
    process_label = option_int(options, "--process_mask_threshold", 0)
    appendix = option_value(options, "--appendix")
    spine_icp_target = option_value(
        options,
        "--spine_icp_target",
        DEFAULT_SPINE_ICP_TARGET,
    )
    spine_reference_path = option_value(
        options,
        "--reference_path",
        str(default_spine_reference_path(spine_icp_target)),
    )
    if args.use_absolute_fe_displacement:
        spine_displacement = absolute_displacement_metadata(
            options,
            str(spine_solve_report_profile(preset=getattr(args, "preset", None))["default_applied_displacement"]),
        )
    else:
        spine_displacement = percent_displacement_metadata(
            model_path,
            report_profile="spine",
            failure_axis="z",
            target_percent=common["solve_and_reporting"]["target_displacement_percent"],
        )
    metadata = common.copy()
    metadata.update({
        "model": "spine-compression",
        "target": {
            "vertebra": appendix,
            "body_label": body_label,
            "process_label": process_label,
        },
        "alignment": {
            "method": SPINE_ALIGNMENT_METHOD,
            "registration_target": spine_icp_target,
            "reference_path": spine_reference_path,
            "registration_scale": option_value(options, "--registration_scale", "auto"),
            "registration_min_scale": option_value(
                options,
                "--registration_min_scale",
                DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
            ),
            "registration_max_scale": option_value(
                options,
                "--registration_max_scale",
                DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
            ),
            "registration_backend": option_value(
                options,
                "--registration_backend",
                DEFAULT_SPINE_REGISTRATION_BACKEND,
            ),
            "registration_landmarks": option_int(
                options,
                "--registration_landmarks",
                DEFAULT_SPINE_REGISTRATION_LANDMARKS,
            ),
            "registration_iterations": option_int(
                options,
                "--registration_iterations",
                DEFAULT_SPINE_REGISTRATION_ITERATIONS,
            ),
        },
        "image_processing": {
            "iso_resolution_mm": option_float(options, "--iso_resolution", DEFAULT_SPINE_ISO_RESOLUTION_MM),
            "registration_label_smoothing": generation_settings["registration_label_smoothing"],
            "spatial_operations": "ICP transform and isotropic output spacing in one shared VTK reslice",
            "image_interpolation": "cubic",
            "label_interpolation": "nearest-neighbor",
            "mask_smoothing": {
                "operation": "binary close/open after ICP resampling",
                "condition": "enabled only when any input spacing dimension exceeds threshold_mm",
                "threshold_mm": option_float(
                    options,
                    "--mask_smoothing_spacing_threshold",
                    DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                ),
            },
            "density_preprocessing": {
                "connectivity_filter": True,
                "bmd_preprocess_threshold": -31,
                "density_binning": {
                    "n_bins": 128,
                    "background_material_id": 0,
                    "trabecular_material_ids": [1, 128],
                    "cortical_material_ids": [129, 256],
                },
            },
        },
        "segmentation": {
            "input_mask": str(args.bone_mask),
            "labels": {
                "vertebral_body": body_label,
                "posterior_process": process_label,
            },
            "derived_masks": {
                "body": "threshold body label",
                "process": "threshold process label",
                "cortical": "density/surface-derived cortical shell generated after alignment",
                "pmma_caps": "fixed-thickness anatomy superior and inferior caps",
            },
        },
        "materials": material_metadata(options, "spine", include_cortical=True),
        "boundary_conditions": {
            "stable_contact": generation_settings["stable_contact"],
            "fixture_geometry": {
                "superior_cap": {
                    "label_id": option_int(options, "--top_node_set_id", 4),
                    "node_set": "body_top",
                    "shape": "fixed-thickness anatomy cap",
                    "pmma_thickness_mm": option_float(
                        options, "--pmma_thick", DEFAULT_SPINE_PMMA_THICKNESS_MM
                    ),
                    "pmma_intrusion_mm": option_float(
                        options, "--pmma_intrusion", DEFAULT_SPINE_PMMA_INTRUSION_MM
                    ),
                    "meaning": (
                        "fixed-thickness anatomy cap: pmma_thickness_mm is total cap thickness; "
                        "pmma_intrusion_mm controls how far anatomy can occupy that fixed "
                        "thickness without overwriting body bone"
                    ),
                },
                "inferior_cap": {
                    "label_id": option_int(options, "--bottom_node_set_id", 3),
                    "node_set": "body_bottom",
                    "shape": "fixed-thickness anatomy cap",
                    "pmma_thickness_mm": option_float(
                        options, "--pmma_thick", DEFAULT_SPINE_PMMA_THICKNESS_MM
                    ),
                    "pmma_intrusion_mm": option_float(
                        options, "--pmma_intrusion", DEFAULT_SPINE_PMMA_INTRUSION_MM
                    ),
                    "meaning": (
                        "fixed-thickness anatomy cap: pmma_thickness_mm is total cap thickness; "
                        "pmma_intrusion_mm controls how far anatomy can occupy that fixed "
                        "thickness without overwriting body bone"
                    ),
                },
            },
            "constraints": [
                {
                    "name": "top_displacement",
                    "node_set": "body_top",
                    "axis": "z",
                    **spine_displacement,
                    "meaning": "superior PMMA cap prescribed toward inferior cap",
                },
                {
                    "name": "bottom_fixed",
                    "node_set": "body_bottom",
                    "axes": ["x", "y", "z"],
                    "value_mm": 0.0,
                    "meaning": "inferior PMMA cap fixed in all displacement directions",
                },
            ],
        },
    })
    return metadata


def _femur_metadata(common, model_path, options, args, generation_settings):
    femur_side = str(option_value(options, "--femur_side", "1"))
    side = "left" if femur_side == "1" else "right" if femur_side == "2" else "unknown"
    compartment_mask = option_path(options, "--compartment_mask")
    if args.use_absolute_fe_displacement:
        femur_displacement = absolute_displacement_metadata(
            options,
            str(femur_solve_report_profile()["default_applied_displacement"]),
        )
    else:
        femur_displacement = percent_displacement_metadata(
            model_path,
            report_profile="femur",
            failure_axis="y",
            target_percent=common["solve_and_reporting"]["target_displacement_percent"],
        )
    metadata = common.copy()
    metadata.update({
        "model": "femur-sideways",
        "target": {
            "side": side,
            "femur_side_code": int(femur_side) if femur_side.isdigit() else femur_side,
            "mask_threshold": option_int(options, "--mask_threshold", 1),
        },
        "alignment": {
            "method": "ICP to side-specific femur reference from cropped/padded input geometry",
            "reference_path": option_value(options, "--reference_path", "bundled side-specific femur reference"),
            "registration_backend": DEFAULT_FEMUR_REGISTRATION_BACKEND,
            "registration_landmarks": DEFAULT_FEMUR_REGISTRATION_LANDMARKS,
            "registration_iterations": DEFAULT_FEMUR_REGISTRATION_ITERATIONS,
        },
        "image_processing": {
            "iso_resolution_mm": option_float(options, "--iso_resolution", DEFAULT_FEMUR_ISO_RESOLUTION_MM),
            "spatial_operations": (
                "ICP transform and isotropic output spacing in one shared VTK reslice"
            ),
            "image_interpolation": "cubic",
            "label_interpolation": "nearest-neighbor",
            "mask_smoothing": {
                "operation": "binary close/open after ICP resampling",
                "condition": "enabled only when any input spacing dimension exceeds threshold_mm",
                "threshold_mm": option_float(
                    options,
                    "--mask_smoothing_spacing_threshold",
                    DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                ),
                "applies_to": [
                    "whole-femur mask",
                    "derived cortical binary mask when compartment_mask is supplied",
                ],
            },
            "density_preprocessing": {
                "connectivity_filter": True,
                "bmd_preprocess_threshold": -31,
                "density_binning": {
                    "n_bins": 128,
                    "background_material_id": 0,
                    "trabecular_material_ids": [1, 128],
                    "cortical_material_ids": [129, 256] if compartment_mask is not None else None,
                },
            },
        },
        "segmentation": {
            "input_mask": str(args.bone_mask),
            "whole_bone_mask_threshold": option_int(options, "--mask_threshold", 1),
            "compartment_mask": compartment_mask,
            "compartment_labels": None
            if compartment_mask is None
            else {
                "cortical": option_int(options, "--cortical_label", 1),
                "trabecular": option_int(options, "--trabecular_label", 2),
            },
        },
        "shaft_standardization": {
            "cut_mode": DEFAULT_FEMUR_CUT_MODE,
            "crop_stage": "fixed rough crop before ICP for registration, final GT-disk-relative flat crop after ICP on the full scan",
            "coverage_definition": "GT disk distal voxel face minus safe flat face clearing the full transformed native scan end",
            "geometry_sidecar": str(model_path.with_name(model_path.stem + "_shaft_geometry.json")),
            "rough_pre_icp_crop": {
                "enabled": True,
                "retained_length_mm": option_float(options, "--femur_shaft_length", 120.0),
                "registration_only": True,
            },
            "greater_trochanter_length": {
                "retained_length_mm": option_float(
                    options,
                    "--femur_greater_trochanter_distal_length",
                    DEFAULT_FEMUR_GREATER_TROCHANTER_DISTAL_LENGTH_MM,
                ),
                "gt_inclusion_length_mm": option_float(
                    options,
                    "--femur_greater_trochanter_inclusion_length",
                    DEFAULT_FEMUR_GREATER_TROCHANTER_INCLUSION_LENGTH_MM,
                ),
                "length_origin": "detected greater-trochanter disk distal edge",
            },
            "cut_plane": "flat post-ICP aligned-frame crop face at fixed shaft length below detected GT-disk distal edge",
            "incomplete_fov_behavior": "fail model generation",
        },
        "materials": material_metadata(options, "femur", include_cortical=compartment_mask is not None),
        "boundary_conditions": {
            "fixture_geometry": {
                "femoral_head": {
                    "node_set": "Femoral_Head_PMMA_Nodes",
                    "shape": "rectangle fixture cap",
                    "relative_to": "model_bbox",
                    "center_fraction": list(FEMORAL_HEAD_FIXTURE_CENTER_FRACTION),
                    "size_fraction": list(SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION),
                    "projection_axis": "y",
                    "pmma_thickness_mm": option_float(options, "--pmma_thick", DEFAULT_PMMA_THICKNESS_MM),
                    "pmma_intrusion_mm": option_float(options, "--pmma_intrusion", DEFAULT_PMMA_INTRUSION_MM),
                    "meaning": "bbox-scaled high-y contact fixture with unsupported columns cropped",
                },
                "greater_trochanter": {
                    "node_set": "Greater_Trochanter_PMMA_Nodes",
                    "shape": "rectangle fixture cap",
                    "relative_to": "model_bbox",
                    "center_fraction": list(GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION),
                    "size_fraction": list(SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION),
                    "projection_axis": "y",
                    "pmma_thickness_mm": option_float(options, "--pmma_thick", DEFAULT_PMMA_THICKNESS_MM),
                    "pmma_intrusion_mm": option_float(options, "--pmma_intrusion", DEFAULT_PMMA_INTRUSION_MM),
                    "meaning": "bbox-scaled low-y contact fixture with unsupported columns cropped",
                },
                "distal_shaft": {
                    "node_set": "Distal_Femur_Nodes",
                    "support_surface": "central 90% straight patch on the post-ICP flat distal shaft crop face",
                    "relative_to": "model_grid",
                    "support_fraction": POST_ICP_DISTAL_SHAFT_SUPPORT_FRACTION,
                    "normal_source": "model-grid distal z direction",
                },
            },
            "constraints": [
                {
                    "name": "top_displacement",
                    "node_set": "Femoral_Head_PMMA_Nodes",
                    "axis": "y",
                    **femur_displacement,
                    "meaning": "femoral head PMMA cap prescribed toward greater trochanter",
                },
                {
                    "name": "bottom_fixed_y_PMMA",
                    "node_set": "Greater_Trochanter_PMMA_Nodes",
                    "axis": "y",
                    "value_mm": 0.0,
                    "meaning": "greater trochanter PMMA cap constrained in loading direction",
                },
                {
                    "name": "bottom_fixed_x",
                    "node_set": "Distal_Femur_Nodes",
                    "axis": "x",
                    "value_mm": 0.0,
                    "meaning": "distal shaft rigid-body constraint",
                },
                {
                    "name": "bottom_fixed_z",
                    "node_set": "Distal_Femur_Nodes",
                    "axis": "z",
                    "value_mm": 0.0,
                    "meaning": "distal shaft rigid-body constraint",
                },
            ],
        },
    })

    return metadata
