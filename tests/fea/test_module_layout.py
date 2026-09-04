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
