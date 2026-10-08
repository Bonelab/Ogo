"""Display smoothing must not alter solved values or blend supports into bone."""

import numpy as np
import pytest


def test_masked_sed_smoothing_excludes_disks_and_preserves_source():
    from ogo.fea.qc_render import smooth_bone_sed

    indices = np.indices((3, 3, 3)).reshape(3, -1).T
    bone = indices[:, 0] < 2
    values = np.where(bone, 4.0, 1000.0)
    original = values.copy()
    display = smooth_bone_sed(values, indices, bone, np.ones(3))
    np.testing.assert_allclose(display, 4)
    np.testing.assert_array_equal(values, original)


def test_sed_missing_or_invalid_is_rejected():
    from ogo.fea.qc_render import smooth_bone_sed

    with pytest.raises(ValueError, match="finite"):
        smooth_bone_sed(np.array([np.nan]), np.zeros((1, 3), dtype=int),
                        np.ones(1, dtype=bool), np.ones(3))


def test_qc_uses_existing_visualizer_extension():
    import inspect
    from ogo.cli.Visualize import vis3d

    signature = inspect.signature(vis3d)
    assert "renderer_setup" in signature.parameters
    assert "label_palette" in signature.parameters


def test_ineligible_qc_is_anatomy_only(tmp_path):
    from ogo.fea.qc_render import export_model_qc

    with pytest.raises(ValueError, match="Ineligible"):
        export_model_qc(tmp_path / "short.n88model", "hip", sed=True, ineligible=True)


def test_render_failure_warns_without_interrupting_solver(monkeypatch):
    from ogo.fea import qc_render

    def fail(*args, **kwargs):
        raise RuntimeError("No graphics display")

    monkeypatch.setattr(qc_render, "export_model_qc", fail)
    with pytest.warns(RuntimeWarning, match="No graphics display"):
        assert qc_render.try_export_model_qc("model.n88model", "hip") is None


def test_too_short_branch_exports_before_exit():
    import inspect
    from ogo.fea import femur

    source = inspect.getsource(femur)
    branch = source.split('if shaft_crop.get("status") == "too_short":', 1)[1].split('sys.exit(1)', 1)[0]
    assert 'model=available_model, ineligible=True' in branch
    assert 'shaft_crop["exports"]' in branch
    assert 'write_model(' not in branch


def test_export_keeps_only_panel_and_settings(tmp_path, monkeypatch):
    import json
    from PIL import Image
    from ogo.fea import qc_render

    def fake_panel(model_path, temporary, site, body_mask, sed, model, ineligible):
        assert ineligible
        assert model is sentinel
        image = temporary / "panel.png"
        Image.new("RGB", (20, 60), "white").save(image)
        for view in ("oblique", "top", "bottom"):
            (temporary / (model_path.stem + "_preview_" + view + ".json")).write_text(
                json.dumps({"eligible_for_solve": False, "view": view}))
        return image

    class Window:
        def GetClassName(self):
            return "vtkEGLRenderWindow"

        def Finalize(self):
            pass

    sentinel = object()
    monkeypatch.setattr(qc_render.vtk, "vtkRenderWindow", Window)
    monkeypatch.setattr(qc_render, "panel", fake_panel)
    path = qc_render.export_model_qc(tmp_path / "short.n88model", "hip",
                                    model=sentinel, ineligible=True)
    assert path.endswith("short_qc_3d.png")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["short_qc_3d.json", "short_qc_3d.png"]
    assert all(not v["eligible_for_solve"] for v in json.loads((tmp_path / "short_qc_3d.json").read_text())["views"])
