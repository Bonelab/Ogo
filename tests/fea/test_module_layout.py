from pathlib import Path


def test_fe_workflow_builders_live_in_domain_modules():
    repo_root = Path(__file__).resolve().parents[2]
    fea_dir = repo_root / "ogo" / "fea"

    assert (fea_dir / "femur.py").exists()
    assert (fea_dir / "spine.py").exists()
    assert not (fea_dir / "sideways_fall.py").exists()
    assert not (fea_dir / "spine_compression.py").exists()


def test_domain_modules_expose_workflow_entry_points():
    from ogo.fea import femur, spine

    assert callable(femur.main)
    assert callable(femur.sidewaysFallFe)
    assert callable(spine.main)
    assert callable(spine.process_vertebra)


def test_femur_exposes_only_current_cropping_helpers():
    from ogo.fea import femur

    assert callable(femur.crop_vtk_images_to_fixed_proximal_length)
    assert callable(femur.crop_vtk_images_to_greater_trochanter_length)
    for name in (
        "crop_vtk_images_to_bbox_ratio",
        "crop_vtk_images_to_proximal_box_ratio",
        "crop_vtk_images_to_flat_post_icp_ratio",
        "crop_vtk_images_to_oblique_post_icp_ratio",
        "standardize_femur_shaft_length",
        "detect_lesser_trochanter_cut_z",
    ):
        assert not hasattr(femur, name), name


def test_fea_shared_dependencies_are_available():
    import inspect
    from ogo.util import Helper
    from ogo.cli.Visualize import vis3d
    from ogo.fea import femur, spine

    for name in ("readNii", "readPolyData", "transformResample",
                 "prepareFiniteElementImage", "imageConnectivity", "bmd_preprocess",
                 "cast2short", "applyMaskByArray", "density2materialID", "maskThreshold"):
        assert callable(getattr(Helper, name)), name
    assert "interpolation" in inspect.signature(Helper.transformResample).parameters
    assert {"renderer_setup", "label_palette"} <= set(inspect.signature(vis3d).parameters)
    data = Path(femur.__file__).parents[1] / "dat"
    assert (data / "LT_FEMUR_SIDEWAYS_FALL_REF.vtk").is_file()
    for target in spine.SPINE_ICP_TARGETS:
        assert spine.default_spine_reference_path(target).is_file()


def test_alignment_and_qc_are_grouped_by_responsibility():
    from ogo.fea.alignment import hip
    from ogo.fea.qc import gallery, images, legacy, render, review, validation

    assert Path(hip.__file__).parent.name == "alignment"
    for module in (gallery, images, legacy, render, review, validation):
        assert Path(module.__file__).parent.name == "qc"
