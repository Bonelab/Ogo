"""Compact QC images preserve dimensions and legacy preview discovery."""

from pathlib import Path

import pytest
from PIL import Image


def test_review_image_is_webp_without_resizing(tmp_path):
    from ogo.fea.qc_images import save_review_image

    path = tmp_path / 'case_qc_3d.webp'
    save_review_image(Image.new('RGB', (1400, 2550), 'white'), path)
    with Image.open(path) as image:
        assert image.format == 'WEBP'
        assert image.size == (1400, 2550)


@pytest.mark.parametrize('extension', ['.png', '.webp'])
def test_gallery_splits_both_formats_into_webp(tmp_path, extension):
    from ogo.fea.gallery import _camera_images

    source = tmp_path / ('case_qc_3d' + extension)
    Image.new('RGB', (1400, 2550), 'white').save(source)
    assets = _camera_images(source, tmp_path / 'gallery.html')
    assert [view for view, _ in assets] == ['oblique', 'top', 'bottom']
    for _, path in assets:
        assert path.suffix == '.webp'
        with Image.open(path) as image:
            assert image.size == (1400, 850)


def test_measurements_discover_legacy_png_and_prefer_webp(tmp_path):
    import csv
    from ogo.fea.validation import write_measurements

    model = tmp_path / 'case.n88model'
    png = tmp_path / 'case_qc_3d.png'
    png.touch()
    path = write_measurements(model, {'site': 'hip'})
    with path.open() as stream:
        assert next(csv.DictReader(stream))['anatomy_image'] == str(png)
    webp = png.with_suffix('.webp')
    webp.touch()
    path = write_measurements(model, {'site': 'hip'})
    with path.open() as stream:
        row = next(csv.DictReader(stream))
    assert row['anatomy_image'] == str(webp)
    assert Path(row['sed_image']).suffix == '.webp'
