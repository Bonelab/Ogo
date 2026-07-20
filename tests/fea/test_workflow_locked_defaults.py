import json
from pathlib import Path

from ogo.fea import femur
from ogo.fea import spine


def test_maintained_fe_workflow_defaults_match_locked_fixture():
    expected_path = Path(__file__).parent / "fixtures" / "fe_workflow_defaults_locked.json"
    expected = json.loads(expected_path.read_text(encoding="utf-8"))

    assert _locked_workflow_defaults() == expected


def _locked_workflow_defaults():
    return {
        "schema_version": 1,
        "workflows": {
            "spine-compression": {
                "iso_resolution_mm": spine.DEFAULT_SPINE_ISO_RESOLUTION_MM,
                "pmma_thickness_mm": spine.DEFAULT_SPINE_PMMA_THICKNESS_MM,
                "pmma_intrusion_mm": spine.DEFAULT_SPINE_PMMA_INTRUSION_MM,
                "pmma_material_id": spine.DEFAULT_SPINE_PMMA_MATERIAL_ID,
                "pmma_e_mpa": spine.DEFAULT_SPINE_PMMA_E_MPA,
                "target_displacement_percent": spine.DEFAULT_SPINE_TARGET_DISPLACEMENT_PERCENT,
                "label_smoothing_sigma_mm": spine.DEFAULT_SPINE_LABEL_SMOOTHING_SIGMA_MM,
                "label_smoothing_threshold": spine.DEFAULT_SPINE_LABEL_SMOOTHING_THRESHOLD,
                "mask_smoothing_spacing_threshold_mm": spine.DEFAULT_SPINE_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                "top_node_set_id": spine.DEFAULT_SPINE_TOP_NODE_SET_ID,
                "bottom_node_set_id": spine.DEFAULT_SPINE_BOTTOM_NODE_SET_ID,
                "registration_backend": spine.DEFAULT_SPINE_REGISTRATION_BACKEND,
                "registration_min_scale": spine.DEFAULT_SPINE_REGISTRATION_MIN_SCALE,
                "registration_max_scale": spine.DEFAULT_SPINE_REGISTRATION_MAX_SCALE,
                "reference_filename": spine.DEFAULT_SPINE_REFERENCE_FILENAME,
                "contact_size_fraction": list(spine.SPINE_CONTACT_SIZE_FRACTION),
                "superior_contact_center_fraction": list(
                    spine.SPINE_SUPERIOR_CONTACT_CENTER_FRACTION
                ),
                "inferior_contact_center_fraction": list(
                    spine.SPINE_INFERIOR_CONTACT_CENTER_FRACTION
                ),
                "benchmark_linear": _selected_spine_benchmark(
                    spine.benchmark_linear_params()
                ),
                "benchmark_nonlinear": _selected_spine_benchmark(
                    spine.benchmark_nonlinear_params()
                ),
            },
            "hip-sideways-fall": {
                "left_femur_code": femur.LEFT_FEMUR,
                "right_femur_code": femur.RIGHT_FEMUR,
                "iso_resolution_mm": femur.DEFAULT_FEMUR_ISO_RESOLUTION_MM,
                "fe_displacement_mm": femur.DEFAULT_FEMUR_FE_DISPLACEMENT,
                "target_displacement_percent": femur.DEFAULT_FEMUR_TARGET_DISPLACEMENT_PERCENT,
                "mask_smoothing_spacing_threshold_mm": femur.DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                "fixture_shape": femur.SIDEWAYS_FALL_FIXTURE_SHAPE,
                "fixture_size_fraction": list(femur.SIDEWAYS_FALL_FIXTURE_SIZE_FRACTION),
                "femoral_head_center_fraction": list(
                    femur.FEMORAL_HEAD_FIXTURE_CENTER_FRACTION
                ),
                "greater_trochanter_center_fraction": list(
                    femur.GREATER_TROCHANTER_FIXTURE_CENTER_FRACTION
                ),
                "distal_shaft_center_fraction": list(
                    femur.DISTAL_SHAFT_FIXTURE_CENTER_FRACTION
                ),
                "distal_shaft_size_fraction": list(
                    femur.DISTAL_SHAFT_FIXTURE_SIZE_FRACTION
                ),
                "distal_shaft_normal": list(femur.DISTAL_SHAFT_FIXTURE_NORMAL),
                "bbox_ratio": list(femur.DEFAULT_FEMUR_BBOX_RATIO),
                "bbox_crop_from": list(femur.DEFAULT_FEMUR_BBOX_CROP_FROM),
                "reference_min_scale": list(femur.DEFAULT_FEMUR_REFERENCE_MIN_SCALE),
                "reference_max_scale": list(femur.DEFAULT_FEMUR_REFERENCE_MAX_SCALE),
                "pmma_thickness_mm": femur.DEFAULT_PMMA_THICKNESS_MM,
                "pmma_intrusion_mm": femur.DEFAULT_PMMA_INTRUSION_MM,
                "input_margin_mm": femur.DEFAULT_FEMUR_INPUT_MARGIN_MM,
                "cut_mode": femur.DEFAULT_FEMUR_CUT_MODE,
                "shaft_length_mm": femur.DEFAULT_FEMUR_SHAFT_LENGTH_MM,
                "lesser_trochanter_distal_offset_mm": femur.DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
                "cortical_label": femur.DEFAULT_CORTICAL_LABEL,
                "trabecular_label": femur.DEFAULT_TRABECULAR_LABEL,
                "node_sets": list(femur.SIDEWAYS_FALL_NODE_SETS),
            },
        },
    }


def _selected_spine_benchmark(params):
    keys = (
        "pmma_mat_id",
        "pmma_E",
        "pmma_v",
        "top_node_set_id",
        "bottom_node_set_id",
        "top_direction",
        "bottom_direction",
        "fe_displacement",
        "target_displacement_percent",
        "elastic_E_func",
        "yield_comp_func",
        "yield_tens_func",
        "cort_elastic_E_func",
        "cort_yield_comp_func",
        "cort_yield_tens_func",
    )
    return {
        key: list(params[key]) if isinstance(params[key], tuple) else params[key]
        for key in keys
    }
