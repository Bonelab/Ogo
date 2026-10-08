from ogo.fea import spine


def test_default_spine_cap_geometry_matches_maintained_settings():
    assert spine.DEFAULT_SPINE_PMMA_THICKNESS_MM == 10
    assert spine.DEFAULT_SPINE_PMMA_INTRUSION_MM == 6
    assert spine.DEFAULT_SPINE_REGISTRATION_BACKEND == "numpy"
    assert spine.DEFAULT_SPINE_REGISTRATION_LANDMARKS == 8000
    assert spine.DEFAULT_SPINE_REGISTRATION_ITERATIONS == 50
    assert spine.DEFAULT_SPINE_ICP_TARGET == "body"


def test_spine_reference_path_depends_on_icp_target():
    assert spine.default_spine_reference_path("body").name == "L4_BODY_SPINE_COMPRESSION_REF.vtk"
    assert spine.default_spine_reference_path("vertebra").name == "L4_FULL_VERTEBRA_SPINE_COMPRESSION_REF.vtk"


def test_benchmark_presets_match_spinefe_notebook_settings():
    linear = spine.benchmark_linear_params()
    nonlinear = spine.benchmark_nonlinear_params()

    assert linear["fe_displacement"] == -0.2
    assert linear["target_displacement_percent"] == 0.68
    assert linear["elastic_E_func"] == "kopperdahl_trab_E"
    assert linear["yield_comp_func"] is None
    assert linear["pmma_yield_compression"] is None
    assert nonlinear["fe_displacement"] == -2.0
    assert nonlinear["target_displacement_percent"] == 4.0
    assert nonlinear["yield_comp_func"] == "kopperdahl_trab_yc"
    assert nonlinear["pmma_yield_compression"] == 70.0
