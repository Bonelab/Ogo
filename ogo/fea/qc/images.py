"""Compact review images; numerical models and source images are not modified."""

from pathlib import Path

REVIEW_PANEL_SIZE = (600, 1092)


def save_qc_panel(image, path):
    """Save the three-view debug panel at review resolution, not publication size."""
    from PIL import Image

    if image.size in ((1400, 2550), REVIEW_PANEL_SIZE):
        image = image.resize(REVIEW_PANEL_SIZE, Image.Resampling.LANCZOS)
    else:
        image = image.copy()
        image.thumbnail(REVIEW_PANEL_SIZE, Image.Resampling.LANCZOS)
    image.convert('RGB').save(path, format='WEBP', quality=85, method=2)


def save_review_image(image, path):
    """Save an RGB WebP at quality 85 without resizing the supplied image."""
    image.convert('RGB').save(path, format='WEBP', quality=85, method=4)


def qc_image_path(model_path, sed=False):
    """Prefer current WebP previews, accepting existing legacy PNGs."""
    model_path = Path(model_path)
    suffix = '_sed_3d' if sed else '_qc_3d'
    path = model_path.with_name(model_path.stem + suffix + '.webp')
    legacy = path.with_suffix('.png')
    return legacy if not path.exists() and legacy.exists() else path
