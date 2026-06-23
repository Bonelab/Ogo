import pytest

from ogo.fea import femur


def test_side_suffix_and_rotation_match_legacy_outputs():
    assert femur.side_suffix(1) == "_LT_FEMUR_SF"
    assert femur.side_suffix(2) == "_RT_FEMUR_SF"
    assert femur.side_rotation(1) == 90
    assert femur.side_rotation(2) == -90


def test_sideways_fall_output_name_preserves_legacy_pattern():
    assert femur.sideways_fall_output_name("density.n88model", 1) == "density_LT_FEMUR_SF.n88model"
    assert femur.sideways_fall_output_name("density.n88model", 2) == "density_RT_FEMUR_SF.n88model"


def test_invalid_femur_side_is_rejected():
    with pytest.raises(ValueError):
        femur.side_suffix(3)


def test_tilted_side_support_vector_has_unit_length():
    np = pytest.importorskip("numpy")

    vector = femur.tilted_side_support_vector(-20)

    assert np.linalg.norm(vector) == pytest.approx(1.0)

