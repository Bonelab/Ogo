# /------------------------------------------------------------------------------+
# | 03-Mar-2026                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
"""Read segmentation labels and expose them through ``ogo check-labels``."""

import argparse
import logging
import os

import numpy as np
import SimpleITK as sitk

import ogo.dat.OgoMasterLabels as lb


logger = logging.getLogger(__name__)


def read_labels(input_image):
    """Return sorted label IDs, including background, after UInt8 reading.

    Accept a filename or path-like object. Preserve the historical UInt8
    conversion, including its limitations for values outside 0–255. Missing
    files raise FileNotFoundError; SimpleITK read errors propagate to the caller.
    This operation does not print or terminate the process.
    """
    path = os.fspath(input_image)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input image is not a file: {path}")
    image = sitk.ReadImage(path, sitk.sitkUInt8)
    array = sitk.GetArrayViewFromImage(image)
    return sorted(np.unique(array).tolist())


def find_missing_labels(present, required):
    """Return absent required IDs, preserving their order and duplicates."""
    return [label for label in required if label not in present]


def add_label_descriptions(labels):
    labels_with_desc = []
    for lab in labels:
        desc = lb.labels_dict.get(lab, {}).get('LABEL', 'unknown label')
        labels_with_desc.append(f"{lab} ({desc})")
    return labels_with_desc


def add_arguments(parser):
    """Define the segmentation input and optional required-label check."""
    parser.add_argument("input_image", help="Input image file (*.nii, *.nii.gz)")
    parser.add_argument("--require", type=int, nargs="*", default=[],
                        help="Required labels (e.g. --require 2 10)")
    parser.add_argument("--include-description", action="store_true",
                        help="Include label descriptions in the output")


def main(argv=None):
    """Report labels and return 0, 1 for a read failure, or 2 for missing IDs."""
    parser = argparse.ArgumentParser(
        prog="ogo check-labels",
        description=(
            "Print labels present in a NIfTI segmentation or check required IDs. "
            "Images are read as UInt8; label values outside 0–255 are converted."
        ),
        epilog=(
            "Examples:\n"
            "  ogo check-labels segmentation.nii.gz\n"
            "  ogo check-labels segmentation.nii.gz --require 2 10\n"
            "An empty --require prints the inventory. Read failures return 1; "
            "missing required labels return 2."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_arguments(parser)
    args = parser.parse_args(argv)

    if not args.input_image.lower().endswith((".nii", ".nii.gz")):
        logger.warning(
            "Input image does not end with .nii or .nii.gz; attempting to read: %s",
            args.input_image,
        )
    logger.info("Reading input image: %s", args.input_image)
    try:
        present = read_labels(args.input_image)
    except FileNotFoundError as error:
        logger.error("%s", error)
        logger.debug("Input validation failed", exc_info=True)
        return 1
    except (OSError, RuntimeError):
        logger.error(
            "Cannot read input image: %s. "
            "Check that it is readable and contains a valid image.",
            args.input_image,
        )
        logger.debug("Image read failed", exc_info=True)
        return 1
    logger.debug("Labels present after UInt8 reading: %s", present)

    if args.require:
        missing = find_missing_labels(present, args.require)
        if missing:
            report = (
                ", ".join(add_label_descriptions(missing))
                if args.include_description else missing
            )
            logger.error("Missing required labels: %s", report)
            return 2
        print("All required labels are present.")
    else:
        report = (
            ", ".join(add_label_descriptions(present))
            if args.include_description else present
        )
        print(f"Labels found: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
