"""Generate and solve finite-element models for spine vertebrae and hip femurs.

This module is intentionally a thin orchestration layer. The scientific model
generation lives in ``ogo.fea``; this command provides the maintained
user-facing entry point for running one or many targets.
"""

import argparse
from collections import namedtuple
import importlib.util
import json
from pathlib import Path
import sys
from typing import Callable, List, Optional, Sequence

from ogo.fea.spine import (
    BENCHMARK_LINEAR_FE_DISPLACEMENT_MM,
    BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM,
    DEFAULT_SPINE_FE_DISPLACEMENT_MM,
    DEFAULT_SPINE_ISO_RESOLUTION_MM,
    DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    DEFAULT_SPINE_PMMA_INTRUSION_MM,
    DEFAULT_SPINE_PMMA_THICKNESS_MM,
    DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
    DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
    DEFAULT_SPINE_REGISTRATION_BACKEND,
    DEFAULT_SPINE_TARGET_DISPLACEMENT_PERCENT,
    SPINE_ALIGNMENT_METHOD,
    default_spine_reference_path,
)
from ogo.fea.femur import (
    DEFAULT_FEMUR_BBOX_CROP_FROM,
    DEFAULT_FEMUR_BBOX_RATIO,
    DEFAULT_FEMUR_CUT_MODE,
    DEFAULT_FEMUR_FE_DISPLACEMENT,
    DEFAULT_FEMUR_ISO_RESOLUTION_MM,
    DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT,
    DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
    DEFAULT_PMMA_INTRUSION_MM,
    DEFAULT_PMMA_THICKNESS_MM,
    DISTAL_SHAFT_FIXTURE_CENTER_FRACTION,
    DISTAL_SHAFT_FIXTURE_SIZE_FRACTION,
    FEMORAL_HEAD_FIXTURE_CENTER_FRACTION,
    GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION,
    SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION,
)


SpineTarget = namedtuple("SpineTarget", ["level", "body_label", "process_label"])

SPINE_PRESETS = {
    "none": [],
    "benchmark-linear": [
        "--fe_displacement",
        str(BENCHMARK_LINEAR_FE_DISPLACEMENT_MM),
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
    ],
    "benchmark-nonlinear": [
        "--fe_displacement",
        str(BENCHMARK_NONLINEAR_FE_DISPLACEMENT_MM),
        "--pmma_yield_compression",
        "70.0",
        "--pmma_yield_tension",
        "70.0",
        "--elastic_E_func",
        "kopperdahl_trab_E",
        "--yield_comp_func",
        "kopperdahl_trab_yc",
        "--yield_tens_func",
        "kopperdahl_trab_yc",
        "--cort_elastic_E_func",
        "kopperdahl_trab_E",
        "--cort_yield_comp_func",
        "kopperdahl_trab_yc",
        "--cort_yield_tens_func",
        "kopperdahl_trab_yc",
    ],
}


def remove_extension(path: Path) -> str:
    """Remove common medical-image extensions while keeping the useful stem."""
    name = path.name
    if name.endswith(".nii.gz"):
        return name[:-7]
    return path.stem


def parse_spine_target(value: str) -> SpineTarget:
    """Parse ``LEVEL:BODY_LABEL:PROCESS_LABEL`` target syntax."""
    parts = [part.strip() for part in value.replace(",", ":").split(":")]
    if len(parts) != 3 or not all(parts):
        raise argparse.ArgumentTypeError(
            "spine targets must use LEVEL:BODY_LABEL:PROCESS_LABEL, for example L1:2:1"
        )
    level, body_label, process_label = parts
    try:
        return SpineTarget(level=level, body_label=int(body_label), process_label=int(process_label))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("body and process labels must be integers") from exc


def expand_sides(sides: Sequence[str]) -> List[str]:
    """Expand the friendly hip side syntax into concrete left/right jobs."""
    expanded = []  # type: List[str]
    for side in sides:
        if side == "both":
            expanded.extend(["left", "right"])
        elif side not in expanded:
            expanded.append(side)
    return expanded


def spine_preset_args(name: str) -> List[str]:
    """Return a copy of the lower-level arguments for one spine preset."""
    return list(SPINE_PRESETS[name])


def build_spine_command(
    *,
    calibrated_image: Path,
    bone_mask: Path,
    target: SpineTarget,
    output_path: Optional[Path],
    extra_args: Sequence[str],
) -> List[str]:
    """Build the lower-level spine FE generator command for one vertebra."""
    cmd = [
        str(calibrated_image),
        str(bone_mask),
        "--mask_threshold",
        str(target.body_label),
        "--process_mask_threshold",
        str(target.process_label),
        "--appendix",
        target.level,
    ]
    if output_path is not None:
        cmd.extend(["--output_path", str(output_path)])
    cmd.extend(extra_args)
    return cmd


def build_femur_command(
    *,
    calibrated_image: Path,
    bone_mask: Path,
    side: str,
    output_path: Optional[Path],
    extra_args: Sequence[str],
) -> List[str]:
    """Build the lower-level femur FE generator command for one side."""
    femur_side = {"left": "1", "right": "2"}[side]
    cmd = [
        str(calibrated_image),
        str(bone_mask),
        "--femur_side",
        femur_side,
    ]
    if output_path is not None:
        cmd.extend(["--output_path", str(output_path)])
    cmd.extend(extra_args)
    return cmd


def expected_spine_model_path(calibrated_image: Path, output_path: Optional[Path], target: SpineTarget) -> Path:
    """Predict the spine model path written by the spine builder."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    return output_dir / "{}_{}.n88model".format(remove_extension(calibrated_image), target.level)


def expected_femur_model_path(calibrated_image: Path, output_path: Optional[Path], side: str) -> Path:
    """Predict the femur model path written by the femur builder."""
    output_dir = output_path if output_path is not None else calibrated_image.parent
    stem = "LF" if side == "left" else "RF"
    return output_dir / "{}_{}.n88model".format(remove_extension(calibrated_image), stem)


def ensure_output_directory(output_path: Optional[Path]) -> None:
    """Create an explicit model output directory before generation starts."""
    if output_path is not None:
        output_path.mkdir(parents=True, exist_ok=True)


def option_value(argv: Sequence[str], option: str, default: Optional[str] = None) -> Optional[str]:
    """Return the final value for an option from an argv-style list."""
    value = default
    for index, token in enumerate(argv):
        if token == option and index + 1 < len(argv):
            value = argv[index + 1]
    return value


def option_float(argv: Sequence[str], option: str, default: float) -> float:
    """Return the final numeric value for an argv-style option."""
    return float(option_value(argv, option, str(default)))


def option_optional_float(argv: Sequence[str], option: str) -> Optional[float]:
    """Return the final numeric value for an optional argv-style option."""
    value = option_value(argv, option)
    return None if value is None else float(value)


def option_int(argv: Sequence[str], option: str, default: int) -> int:
    """Return the final integer value for an argv-style option."""
    return int(option_value(argv, option, str(default)))


def option_present(argv: Sequence[str], option: str) -> bool:
    """Return whether a flag-like option is present."""
    return option in argv


def option_path(argv: Sequence[str], option: str, default: Optional[Path] = None) -> Optional[str]:
    """Return the final path value for an option as a string."""
    value = option_value(argv, option, str(default) if default is not None else None)
    return None if value is None else str(value)


def option_n_values(argv: Sequence[str], option: str, count: int, default: Sequence) -> list:
    """Return a fixed-width option value list from an argv-style list."""
    values = list(default)
    for index, token in enumerate(argv):
        if token == option and index + count < len(argv):
            values = list(argv[index + 1 : index + 1 + count])
    parsed = []
    for value in values:
        token = str(value).strip().lower() if value is not None else "none"
        parsed.append(None if token in {"", "none", "null", "auto"} else value)
    return parsed


def target_displacement_percent(args: argparse.Namespace, model_type: str) -> float:
    """Return the maintained endpoint displacement percentage for this model."""
    if args.target_displacement is not None:
        return args.target_displacement
    if model_type == "spine":
        return DEFAULT_SPINE_TARGET_DISPLACEMENT_PERCENT
    return DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT


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


def solve_model(
    model_path: Path,
    args: argparse.Namespace,
    model_type: str,
    generator_argv: Sequence[str],
) -> None:
    """Run the generic FAIM adapter on one generated model."""
    try:
        from ogo.util.faim import run_faim_pipeline
    except ModuleNotFoundError:
        # Support direct script use: python ogo/cli/GenerateFEM.py ...
        module_path = Path(__file__).resolve().parents[1] / "util" / "faim.py"
        spec = importlib.util.spec_from_file_location("ogo_util_faim", str(module_path))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        run_faim_pipeline = module.run_faim_pipeline

    # Spine compression uses z reaction force; sideways-fall hip uses y.
    if model_type == "spine":
        analysis_var = "fz_ns1"
        pistoia_vars = []
        failure_axis = "z"
    else:
        analysis_var = "fy_ns1"
        pistoia_vars = ["pis_fy_fail", "pis_stiffy"]
        failure_axis = "y"

    default_applied_displacement = (
        str(DEFAULT_SPINE_FE_DISPLACEMENT_MM)
        if model_type == "spine"
        else str(DEFAULT_FEMUR_FE_DISPLACEMENT)
    )
    applied_displacement = option_value(
        generator_argv,
        "--fe_displacement",
        default_applied_displacement,
    )
    target_displacement = target_displacement_percent(args, model_type)
    report_profile = "spine" if model_type == "spine" else "femur"

    run_faim_pipeline(
        model_file=model_path,
        output_prefix=model_path.with_suffix(""),
        analysis_var=analysis_var,
        pistoia_vars=pistoia_vars,
        failure_axis=failure_axis,
        threads=args.threads,
        conda_env=args.faim_env,
        conda_executable=args.conda_executable,
        install_root=args.faim_install_root,
        bin_dir=args.faim_bin_dir,
        license_dir=args.faim_license_dir,
        faim_command=args.faim_command,
        n88modelinfo_command=args.n88modelinfo_command,
        n88derivedfields_command=args.n88derivedfields_command,
        n88postfaim_command=args.n88postfaim_command,
        n88pistoia_command=args.n88pistoia_command,
        n88tabulate_command=args.n88tabulate_command,
        n88copymodel_command=args.n88copymodel_command,
        critical_volume=args.critical_volume,
        critical_strain=args.critical_strain,
        exclude=args.exclude,
        run_pistoia=False,
        applied_displacement=applied_displacement,
        target_displacement=target_displacement,
        report_profile=report_profile,
        solve_displacement_percent=target_displacement,
        compress=not args.no_compress,
        require_pistoia=args.require_pistoia,
        dry_run=args.dry_run,
    )


def write_modeling_metadata(
    model_path: Path,
    model_type: str,
    generator_argv: Sequence[str],
    args: argparse.Namespace,
    bc_audit_summary: Optional[dict] = None,
) -> Optional[Path]:
    """Write a traceable model-building record next to a generated n88model."""
    if args.dry_run:
        return None
    if not model_path.exists():
        return None

    output_path = model_path.with_name(model_path.with_suffix("").name + "_modeling.json")
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
            "target_displacement_percent": target_displacement_percent(args, model_type),
            "target_displacement_definition": "percent strain converted from model geometry for spine; percent of femur length for femur",
            "run_pistoia": False,
            "compress_solved_model": not args.no_compress,
        },
    }

    if model_type == "spine":
        body_label = option_int(generator_argv, "--mask_threshold", 0)
        process_label = option_int(generator_argv, "--process_mask_threshold", 0)
        appendix = option_value(generator_argv, "--appendix")
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
                "reference_path": option_value(
                    generator_argv,
                    "--reference_path",
                    str(default_spine_reference_path()),
                ),
                "registration_scale": option_value(generator_argv, "--registration_scale", "auto"),
                "registration_min_scale": option_value(
                    generator_argv,
                    "--registration_min_scale",
                    DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
                ),
                "registration_max_scale": option_value(
                    generator_argv,
                    "--registration_max_scale",
                    DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
                ),
                "registration_backend": option_value(
                    generator_argv,
                    "--registration_backend",
                    DEFAULT_SPINE_REGISTRATION_BACKEND,
                ),
            },
            "image_processing": {
                "iso_resolution_mm": option_float(generator_argv, "--iso_resolution", DEFAULT_SPINE_ISO_RESOLUTION_MM),
                "spatial_operations": "ICP transform and isotropic output spacing in one shared VTK reslice",
                "image_interpolation": "cubic",
                "label_interpolation": "nearest-neighbor",
                "mask_smoothing": {
                    "operation": "binary close/open after ICP resampling",
                    "condition": "enabled only when any input spacing dimension exceeds threshold_mm",
                    "threshold_mm": option_float(
                        generator_argv,
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
            "materials": {
                "builder": "ogo.fea.materials.build_spine_material_table",
                "convention": "region material IDs; 0 background, 1..128 trabecular, 129..256 cortical, PMMA explicit",
                "poissons_ratio": option_float(generator_argv, "--poissons_ratio", 0.3),
                "trabecular": {
                    "elastic_E_func": option_value(generator_argv, "--elastic_E_func", "default_E"),
                    "yield_comp_func": option_value(generator_argv, "--yield_comp_func"),
                    "yield_tens_func": option_value(generator_argv, "--yield_tens_func"),
                    "material_id_range": [1, 128],
                },
                "cortical": {
                    "elastic_E_func": option_value(
                        generator_argv,
                        "--cort_elastic_E_func",
                        option_value(generator_argv, "--elastic_E_func", "default_E"),
                    ),
                    "yield_comp_func": option_value(
                        generator_argv,
                        "--cort_yield_comp_func",
                        option_value(generator_argv, "--yield_comp_func"),
                    ),
                    "yield_tens_func": option_value(
                        generator_argv,
                        "--cort_yield_tens_func",
                        option_value(generator_argv, "--yield_tens_func"),
                    ),
                    "poissons_ratio": option_float(
                        generator_argv,
                        "--cort_poissons_ratio",
                        option_float(generator_argv, "--poissons_ratio", 0.3),
                    ),
                    "material_id_range": [129, 256],
                },
                "pmma": {
                    "material_id": option_int(generator_argv, "--pmma_mat_id", 5000),
                    "elastic_E_MPa": option_float(generator_argv, "--pmma_E", 2500),
                    "poissons_ratio": option_float(generator_argv, "--pmma_v", 0.3),
                    "yield_compression_MPa": option_optional_float(generator_argv, "--pmma_yield_compression"),
                    "yield_tension_MPa": option_optional_float(generator_argv, "--pmma_yield_tension"),
                },
            },
            "boundary_conditions": {
                "fixture_geometry": {
                    "superior_cap": {
                        "label_id": option_int(generator_argv, "--top_node_set_id", 4),
                        "node_set": "body_top",
                        "shape": "fixed-thickness anatomy cap",
                        "pmma_thickness_mm": option_float(
                            generator_argv, "--pmma_thick", DEFAULT_SPINE_PMMA_THICKNESS_MM
                        ),
                        "pmma_intrusion_mm": option_float(
                            generator_argv, "--pmma_intrusion", DEFAULT_SPINE_PMMA_INTRUSION_MM
                        ),
                        "meaning": (
                            "fixed-thickness anatomy cap: pmma_thickness_mm is total cap thickness; "
                            "pmma_intrusion_mm controls how far anatomy can occupy that fixed "
                            "thickness without overwriting body bone"
                        ),
                    },
                    "inferior_cap": {
                        "label_id": option_int(generator_argv, "--bottom_node_set_id", 3),
                        "node_set": "body_bottom",
                        "shape": "fixed-thickness anatomy cap",
                        "pmma_thickness_mm": option_float(
                            generator_argv, "--pmma_thick", DEFAULT_SPINE_PMMA_THICKNESS_MM
                        ),
                        "pmma_intrusion_mm": option_float(
                            generator_argv, "--pmma_intrusion", DEFAULT_SPINE_PMMA_INTRUSION_MM
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
    else:
        femur_side = option_value(generator_argv, "--femur_side", "1")
        side = "left" if femur_side == "1" else "right" if femur_side == "2" else "unknown"
        compartment_mask = option_path(generator_argv, "--compartment_mask")
        femur_cut_mode = option_value(generator_argv, "--femur_cut_mode", DEFAULT_FEMUR_CUT_MODE)
        femur_bbox_ratio = option_n_values(generator_argv, "--femur_bbox_ratio", 3, DEFAULT_FEMUR_BBOX_RATIO)
        femur_bbox_crop_from = option_n_values(
            generator_argv,
            "--femur_bbox_crop_from",
            3,
            DEFAULT_FEMUR_BBOX_CROP_FROM,
        )
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
                "mask_threshold": option_int(generator_argv, "--mask_threshold", 1),
            },
            "alignment": {
                "method": "ICP to side-specific femur reference from cropped/padded input geometry",
                "reference_path": option_value(generator_argv, "--reference_path", "bundled side-specific femur reference"),
            },
            "image_processing": {
                "iso_resolution_mm": option_float(generator_argv, "--iso_resolution", DEFAULT_FEMUR_ISO_RESOLUTION_MM),
                "spatial_operations": (
                    "ICP transform and isotropic output spacing in one shared VTK reslice"
                ),
                "image_interpolation": "cubic",
                "label_interpolation": "nearest-neighbor",
                "mask_smoothing": {
                    "operation": "binary close/open after ICP resampling",
                    "condition": "enabled only when any input spacing dimension exceeds threshold_mm",
                    "threshold_mm": option_float(
                        generator_argv,
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
                "whole_bone_mask_threshold": option_int(generator_argv, "--mask_threshold", 1),
                "compartment_mask": compartment_mask,
                "compartment_labels": None
                if compartment_mask is None
                else {
                    "cortical": option_int(generator_argv, "--cortical_label", 1),
                    "trabecular": option_int(generator_argv, "--trabecular_label", 2),
                },
            },
            "shaft_standardization": {
                "cut_mode": femur_cut_mode,
                "fixed_length_mm": option_float(generator_argv, "--femur_shaft_length", 100.0),
                "lesser_trochanter_distal_offset_mm": option_float(
                    generator_argv,
                    "--femur_lesser_trochanter_distal_offset",
                    DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
                ),
                "bbox_ratio": femur_bbox_ratio,
                "bbox_crop_from": femur_bbox_crop_from,
                "cut_plane": (
                    "input bbox-ratio crop face transformed with the femur"
                    if femur_cut_mode == "bbox_ratio"
                    else "flat model-grid z plane"
                ),
                "incomplete_fov_behavior": "fail model generation",
            },
            "materials": {
                "builder": "ogo.fea.materials.build_femur_material_table",
                "convention": "region material IDs; 0 background, 1..128 trabecular, optional 129..256 cortical, PMMA explicit",
                "include_cortical_region": compartment_mask is not None,
                "poissons_ratio": option_float(generator_argv, "--poissons_ratio", 0.3),
                "trabecular": {
                    "elastic_E_func": option_value(generator_argv, "--elastic_E_func", "default_E"),
                    "yield_comp_func": option_value(generator_argv, "--yield_comp_func"),
                    "yield_tens_func": option_value(generator_argv, "--yield_tens_func"),
                    "material_id_range": [1, 128],
                },
                "cortical": None
                if compartment_mask is None
                else {
                    "elastic_E_func": option_value(
                        generator_argv,
                        "--cort_elastic_E_func",
                        option_value(generator_argv, "--elastic_E_func", "default_E"),
                    ),
                    "yield_comp_func": option_value(
                        generator_argv,
                        "--cort_yield_comp_func",
                        option_value(generator_argv, "--yield_comp_func"),
                    ),
                    "yield_tens_func": option_value(
                        generator_argv,
                        "--cort_yield_tens_func",
                        option_value(generator_argv, "--yield_tens_func"),
                    ),
                    "poissons_ratio": option_float(
                        generator_argv,
                        "--cort_poissons_ratio",
                        option_float(generator_argv, "--poissons_ratio", 0.3),
                    ),
                    "material_id_range": [129, 256],
                },
                "pmma": {
                    "material_id": option_int(generator_argv, "--pmma_mat_id", 5000),
                    "elastic_E_MPa": option_float(generator_argv, "--pmma_E", 2500),
                    "poissons_ratio": option_float(generator_argv, "--pmma_v", 0.3),
                    "yield_compression_MPa": option_optional_float(generator_argv, "--pmma_yield_compression"),
                    "yield_tension_MPa": option_optional_float(generator_argv, "--pmma_yield_tension"),
                },
            },
            "boundary_conditions": {
                "fixture_geometry": {
                    "femoral_head": {
                        "node_set": "Femoral_Head_PMMA_Nodes",
                        "shape": "rectangle fixture cap",
                        "relative_to": "model_bbox",
                        "center_fraction": list(FEMORAL_HEAD_FIXTURE_CENTER_FRACTION),
                        "size_fraction": list(SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION),
                        "projection_axis": "y",
                        "pmma_thickness_mm": option_float(generator_argv, "--pmma_thick", DEFAULT_PMMA_THICKNESS_MM),
                        "pmma_intrusion_mm": option_float(generator_argv, "--pmma_intrusion", DEFAULT_PMMA_INTRUSION_MM),
                        "meaning": "bbox-scaled high-y contact fixture with unsupported columns cropped",
                    },
                    "greater_trochanter": {
                        "node_set": "Greater_Trochanter_PMMA_Nodes",
                        "shape": "rectangle fixture cap",
                        "relative_to": "model_bbox",
                        "center_fraction": list(GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION),
                        "size_fraction": list(SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION),
                        "projection_axis": "y",
                        "pmma_thickness_mm": option_float(generator_argv, "--pmma_thick", DEFAULT_PMMA_THICKNESS_MM),
                        "pmma_intrusion_mm": option_float(generator_argv, "--pmma_intrusion", DEFAULT_PMMA_INTRUSION_MM),
                        "meaning": "bbox-scaled low-y contact fixture with unsupported columns cropped",
                    },
                    "distal_shaft": {
                        "node_set": "Distal_Femur_Nodes",
                        "support_surface": (
                            "finite bbox-relative patch projected onto the transformed oblique shaft surface"
                            if femur_cut_mode == "bbox_ratio"
                            else "flat distal shaft cut face"
                        ),
                        "relative_to": "model_bbox" if femur_cut_mode == "bbox_ratio" else None,
                        "center_fraction": (
                            list(DISTAL_SHAFT_FIXTURE_CENTER_FRACTION)
                            if femur_cut_mode == "bbox_ratio"
                            else None
                        ),
                        "size_fraction": (
                            list(DISTAL_SHAFT_FIXTURE_SIZE_FRACTION)
                            if femur_cut_mode == "bbox_ratio"
                            else None
                        ),
                        "normal_source": (
                            "transformed input bbox-ratio crop face"
                            if femur_cut_mode == "bbox_ratio"
                            else "model-grid distal z face"
                        ),
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

    output_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {output_path}")
    return output_path


def audit_generated_model(model_path: Path, model_type: str, args: argparse.Namespace) -> Optional[dict]:
    """Return BC audit summary and optionally write a debug PNG."""
    if args.skip_bc_audit or args.dry_run:
        return None
    if not model_path.exists():
        return None

    from ogo.cli.CheckFEModelBC import audit_model

    audit_kind = "spine-compression" if model_type == "spine" else "femur-sideways"
    result = audit_model(
        model_path,
        model=audit_kind,
        flat_tolerance=args.bc_audit_flat_tolerance,
        write_json=False,
        write_csv_file=False,
        write_plot=args.debug,
    )
    if result["png_path"] is not None:
        print(f"Wrote {result['png_path']}")
    if not result["passed"]:
        failed = [check["name"] for check in result["checks"] if not check["passed"]]
        raise RuntimeError("BC audit failed for {}: {}".format(model_path, "; ".join(failed)))
    return result["summary"]


def _call_cli(main_func: Callable[[], None], program: str, argv: Sequence[str]) -> None:
    old_argv = sys.argv
    try:
        sys.argv = [program, *argv]
        try:
            main_func()
        except SystemExit as exc:
            if exc.code not in (None, 0):
                raise
    finally:
        sys.argv = old_argv


def run_spine_command(argv: Sequence[str]) -> None:
    from ogo.fea.spine import main as spine_main

    _call_cli(spine_main, "ogoFEA-spine-builder", argv)


def run_femur_command(argv: Sequence[str]) -> None:
    from ogo.fea.femur import main as femur_main

    _call_cli(femur_main, "ogoFEA-hip-builder", argv)


def print_dry_run(program: str, argv: Sequence[str]) -> None:
    print(" ".join([program, *argv]))


def _add_common_image_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("calibrated_image", type=Path, help="Calibrated density image.")
    parser.add_argument("bone_mask", type=Path, help="Bone mask or labelled bone mask.")
    parser.add_argument(
        "--output_path",
        type=Path,
        default=None,
        help="Directory for generated .n88model files. Defaults to the input image directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print generated lower-level commands without running model generation or solving.",
    )
    parser.add_argument(
        "--no-solve",
        action="store_true",
        help="Only generate the .n88model; skip FAIM solve and postprocessing.",
    )
    parser.add_argument("--threads", type=int, default=4, help="FAIM solver thread count.")
    parser.add_argument("--faim_env", default=None, help="Optional conda environment for FAIM/N88 tools.")
    parser.add_argument("--conda_executable", default="conda", help="Conda executable for --faim_env.")
    parser.add_argument("--faim_install_root", default=None, help="Optional FAIM install root.")
    parser.add_argument("--faim_bin_dir", default=None, help="Optional directory containing FAIM/N88 commands.")
    parser.add_argument("--faim_license_dir", default=None, help="Optional Numerics88 license directory.")
    parser.add_argument("--faim_command", default=None, help="Override FAIM solver command.")
    parser.add_argument("--n88modelinfo_command", default=None, help="Override n88modelinfo command.")
    parser.add_argument("--n88derivedfields_command", default=None, help="Override n88derivedfields command.")
    parser.add_argument("--n88postfaim_command", default=None, help="Override n88postfaim command.")
    parser.add_argument("--n88pistoia_command", default=None, help="Override n88pistoia command.")
    parser.add_argument("--n88tabulate_command", default=None, help="Override n88tabulate command.")
    parser.add_argument("--n88copymodel_command", default=None, help="Override n88copymodel command.")
    parser.add_argument(
        "--critical_volume",
        type=float,
        default=2.0,
        help="Pistoia critical volume percentage.",
    )
    parser.add_argument(
        "--critical_strain",
        type=float,
        default=0.007,
        help="Pistoia critical EES strain.",
    )
    parser.add_argument("--exclude", type=int, default=5000, help="Material ID excluded by Pistoia.")
    parser.add_argument(
        "--target_displacement",
        type=float,
        default=None,
        help=(
            "Profile reporting endpoint as percent strain. Defaults to 2.0 for spine "
            "and 4.0 for femur; converted to mm from the generated model geometry."
        ),
    )
    parser.add_argument(
        "--require_pistoia",
        action="store_true",
        help="Fail if Pistoia postprocessing fails.",
    )
    parser.add_argument(
        "--no_compress",
        action="store_true",
        help="Skip n88copymodel --compress after solving.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Write debug sidecars such as BC audit PNGs and spine QC quick looks.",
    )
    parser.add_argument(
        "--skip_bc_audit",
        action="store_true",
        help="Skip automatic boundary-condition audit files after model generation.",
    )
    parser.add_argument(
        "--bc_audit_flat_tolerance",
        type=float,
        default=1.0e-4,
        help="Flat-plane tolerance for automatic boundary-condition audit checks.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ogoFEA",
        description=(
            "Generate N88 finite-element models for spine vertebrae and hip femurs. "
            "Unknown options are forwarded to the selected lower-level FE generator."
        ),
    )
    subparsers = parser.add_subparsers(dest="model_type")

    spine = subparsers.add_parser(
        "spine",
        help="Generate one or more vertebral compression FE models.",
        description=(
            "Generate spine FE models. Repeat --vertebra for each level to run. "
            "Each target uses LEVEL:BODY_LABEL:PROCESS_LABEL."
        ),
    )
    _add_common_image_args(spine)
    spine.add_argument(
        "--vertebra",
        action="append",
        type=parse_spine_target,
        required=True,
        metavar="LEVEL:BODY_LABEL:PROCESS_LABEL",
        help="Vertebral target, for example L1:2:1. Repeat for all levels to process.",
    )
    spine.add_argument(
        "--preset",
        choices=sorted(SPINE_PRESETS),
        default="benchmark-linear",
        help=(
            "Spine FE parameter preset. The default benchmark-linear preset uses the "
            "public spineFE-benchmark linear model settings; use none to pass only "
            "explicit lower-level options."
        ),
    )

    hip = subparsers.add_parser(
        "hip",
        help="Generate hip sideways-fall FE models for left, right, or both femurs.",
    )
    _add_common_image_args(hip)
    hip.add_argument(
        "--side",
        action="append",
        choices=["left", "right", "both"],
        default=None,
        help="Femur side to process. Repeatable; defaults to both.",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = build_parser()
    args, extra_args = parser.parse_known_args(argv)

    if args.model_type is None:
        parser.error("model type is required: choose spine or hip")

    if not args.dry_run:
        ensure_output_directory(args.output_path)

    if args.model_type == "spine":
        spine_extra_args = spine_preset_args(args.preset) + list(extra_args) + [
            "--quality_control",
            str(bool(args.debug)),
        ]
        for target in args.vertebra:
            cmd = build_spine_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                target=target,
                output_path=args.output_path,
                extra_args=spine_extra_args,
            )
            if args.dry_run:
                print_dry_run("ogoFEA-spine-builder", cmd)
            else:
                run_spine_command(cmd)
                model_path = expected_spine_model_path(args.calibrated_image, args.output_path, target)
                bc_audit_summary = audit_generated_model(
                    model_path,
                    "spine",
                    args,
                )
                write_modeling_metadata(model_path, "spine", cmd, args, bc_audit_summary)
            if not args.no_solve:
                solve_model(
                    expected_spine_model_path(args.calibrated_image, args.output_path, target),
                    args,
                    "spine",
                    cmd,
                )
        return

    if args.model_type == "hip":
        for side in expand_sides(args.side or ["both"]):
            cmd = build_femur_command(
                calibrated_image=args.calibrated_image,
                bone_mask=args.bone_mask,
                side=side,
                output_path=args.output_path,
                extra_args=extra_args,
            )
            if args.dry_run:
                print_dry_run("ogoFEA-hip-builder", cmd)
            else:
                run_femur_command(cmd)
                model_path = expected_femur_model_path(args.calibrated_image, args.output_path, side)
                bc_audit_summary = audit_generated_model(
                    model_path,
                    "hip",
                    args,
                )
                write_modeling_metadata(model_path, "hip", cmd, args, bc_audit_summary)
            if not args.no_solve:
                solve_model(
                    expected_femur_model_path(args.calibrated_image, args.output_path, side),
                    args,
                    "hip",
                    cmd,
                )
        return

    parser.error(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    main()
