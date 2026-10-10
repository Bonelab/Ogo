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
