from ogo.fea import spine


def test_default_spine_cap_geometry_matches_maintained_settings():
    assert spine.DEFAULT_SPINE_PMMA_THICKNESS_MM == 10
    assert spine.DEFAULT_SPINE_PMMA_INTRUSION_MM == 6
    assert spine.DEFAULT_SPINE_REGISTRATION_BACKEND == "vtk"


def test_benchmark_presets_match_spinefe_notebook_settings():
    linear = spine.benchmark_linear_params()
    nonlinear = spine.benchmark_nonlinear_params()

    assert linear["fe_displacement"] == -0.2
    assert linear["elastic_E_func"] == "kopperdahl_trab_E"
    assert linear["yield_comp_func"] is None
    assert linear["pmma_yield_compression"] is None
    assert nonlinear["fe_displacement"] == -2.0
    assert nonlinear["yield_comp_func"] == "kopperdahl_trab_yc"
    assert nonlinear["pmma_yield_compression"] == 70.0
