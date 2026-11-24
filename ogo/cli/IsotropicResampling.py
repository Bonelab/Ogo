# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""Resample an image file through ``ogo resample-image``."""

import argparse
import logging
import math
import os

import SimpleITK as sitk

import ogo.util.Helper as ogo


logger = logging.getLogger(__name__)


def resample_image_file(input_filename, output_filename, iso_resolution=1.0,
                        image_type="ct", overwrite=False):
    """Write a resampled image and return its output path as a string.

    Preserve the existing 3D CT/mask calculation and pixel type. Resolution is
    in millimetres. Raise ValueError for invalid settings, FileExistsError for
    an unapproved destination, and file/SimpleITK errors for I/O failures.
    Never prompt or exit. A later write failure does not promise rollback.
    """
    if not math.isfinite(iso_resolution) or iso_resolution <= 0:
        raise ValueError(
            "Isotropic resolution must be finite and greater than zero")
    if image_type not in ("ct", "mask"):
        raise ValueError(f"Unknown image type: {image_type!r}; use ct or mask")
    source, target = os.fspath(input_filename), os.fspath(output_filename)
    if os.path.lexists(target) and not overwrite:
        raise FileExistsError(
            f"Output already exists: {target}. Use --overwrite to replace it")
    if not os.path.isfile(source):
        raise FileNotFoundError(f"Input image is not a file: {source}")
    logger.info("Reading input image: %s", source)
    image = sitk.ReadImage(source)
    output = ogo.isotropicResampling(image, iso_resolution, image_type)
    logger.info("Writing output image: %s", target)
    sitk.WriteImage(output, target)
    logger.info("Completed: %s", target)
    return target


def add_arguments(parser):
    """Define input/output paths and existing resampling options."""
    parser.add_argument("input_filename", help="Input image (*.nii, *.nii.gz)")
    parser.add_argument("output_filename",
                        help="Output image (*.nii, *.nii.gz)")
    parser.add_argument("--iso-resolution", type=float, default=1.0, metavar="MM",
                        help="Isotropic resolution in mm (default: %(default)s)")
    parser.add_argument("--image-type", choices=["ct", "mask"], default="ct",
                        help="Input type (default: %(default)s)")
    parser.add_argument("-o", "--overwrite", action="store_true",
                        help="Permit replacement of an existing output")


def main(argv=None):
    """Parse arguments, resample, and return success or failure."""
    parser = argparse.ArgumentParser(
        prog="ogo resample-image",
        description="Resample a 3D CT or mask to isotropic spacing. CT uses existing "
                    "Gaussian antialiasing and B-spline interpolation; masks use nearest neighbor.",
        epilog="Existing outputs require --overwrite. No keyboard prompts are used.",
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    logger.info("Parameters:\n  Input:      %s\n  Output:     %s\n  Resolution: %s mm\n  Image type: %s\n  Overwrite:  %s",
                args.input_filename, args.output_filename, args.iso_resolution,
                args.image_type, args.overwrite)
    try:
        resample_image_file(**vars(args))
    except (OSError, RuntimeError, ValueError) as error:
        logger.error("Resampling failed: %s", error)
        logger.debug("Resampling failure details", exc_info=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ===============================================================================
# Retained old resample function if needed:
#    # Determine sigma
#    ogo.message('Resampling to a resolution of {}'.format(iso_resolution))
#    sigma = iso_resolution * np.sqrt(2 * np.log(2)) / np.pi
#    ogo.message('  For a half-power cut-off frequency of {:0.6f}, '
#               'a standard deviation of {:0.6f} is being used.'.format(1/(2.0*iso_resolution), sigma))
#
#    # Filter each sampling directions as needed
#    for i, spacing in enumerate(image.GetSpacing()):
#        if spacing > iso_resolution:
#            ogo.message('  No antialiasing in direction {}'.format(i))
#            continue
#
#        ogo.message('  Antialiasing in direction {}'.format(i))
#        image = sitk.RecursiveGaussian(
#          image,
#          sigma, False,
#          sitk.RecursiveGaussianImageFilter.ZeroOrder,
#          i
#        )
#
#    # Determine the output size
#    resolution = [iso_resolution for i in range(image.GetDimension())]
#    size = [int(np.ceil(s * i / o)) for o, i, s in
#            zip(resolution, image.GetSpacing(), image.GetSize())]
#    ogo.message('  Input Size:  {}'.format(image.GetSize()))
#    ogo.message('  Output Size: {}'.format(size))
#
#    ogo.message('Finding minimum in image')
#    stats = sitk.StatisticsImageFilter()
#    stats.Execute(image)
#    vox_min = stats.GetMinimum()
#    ogo.message('  Minimum intensity: {}'.format(vox_min))
#
#    ogo.message('Resampling...')
#    transform = sitk.Euler3DTransform()
#    transform.SetIdentity()
#
#    ogo.message('  Using BSpline interpolation')
#
#    output = sitk.Resample(
#      image,
#      size,
#      transform,
#      sitk.sitkBSpline,
#      image.GetOrigin(),
#      resolution,
#      image.GetDirection(),
#      vox_min,
#      image.GetPixelID()
#    )
