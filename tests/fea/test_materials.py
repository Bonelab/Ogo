import pytest

from ogo.fea import material_laws
from ogo.fea import materials


def test_resolve_material_func_accepts_callable_name_and_none():
    def custom_law(density):
        return density

    assert materials.resolve_material_func(custom_law) is custom_law
    assert materials.resolve_material_func("kopperdahl_trab_E") is material_laws.kopperdahl_trab_E
    assert materials.resolve_material_func(None) is None


def test_spine_material_log_values_defaults_to_default_E():
    values = materials.spine_material_log_values({})

    assert values["elastic_E_func"] is material_laws.default_E
    assert values["cort_elastic_E_func"] is material_laws.default_E
    assert values["yield_comp_func"] is None
    assert values["cort_poissons_ratio"] == 0.3


def test_spine_material_log_values_resolves_benchmark_laws():
    values = materials.spine_material_log_values(
        {
            "elastic_E_func": "kopperdahl_trab_E",
            "yield_comp_func": "kopperdahl_trab_yc",
            "cort_elastic_E_func": "kopperdahl_trab_E",
        }
    )

    assert values["elastic_E_func"] is material_laws.kopperdahl_trab_E
    assert values["yield_comp_func"] is material_laws.kopperdahl_trab_yc
    assert values["cort_elastic_E_func"] is material_laws.kopperdahl_trab_E


def test_build_spine_material_table_requires_vtkbone_when_exercised():
    pytest.importorskip("vtkbone")

    table = materials.build_spine_material_table(
        [100.0, 200.0],
        n_bins=2,
        pmma_mat_id=300,
        elastic_E_func="kopperdahl_trab_E",
        cort_elastic_E_func="kopperdahl_trab_E",
    )

    assert table.GetNumberOfMaterials() == 5


def test_build_femur_material_table_can_use_single_trabecular_region():
    pytest.importorskip("vtkbone")

    table = materials.build_femur_material_table(
        [100.0, 200.0],
        n_bins=2,
        pmma_mat_id=300,
        elastic_E_func="kopperdahl_trab_E",
    )

    assert table.GetNumberOfMaterials() == 3


def test_build_femur_material_table_can_use_trab_cort_regions():
    pytest.importorskip("vtkbone")

    table = materials.build_femur_material_table(
        [100.0, 200.0, 100.0, 200.0],
        n_bins=2,
        pmma_mat_id=300,
        elastic_E_func="kopperdahl_trab_E",
        cort_elastic_E_func="bayraktar_cort_E",
        include_cortical=True,
    )

    assert table.GetNumberOfMaterials() == 5
