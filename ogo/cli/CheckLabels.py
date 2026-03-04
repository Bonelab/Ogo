# /------------------------------------------------------------------------------+
# | 03-Mar-2026                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
import argparse
import os
import sys
import SimpleITK as sitk
import numpy as np
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def add_label_descriptions(labels):
    labels_with_desc = []
    for lab in labels:
        desc = lb.labels_dict.get(lab, {}).get('LABEL', 'unknown label')
        labels_with_desc.append(f"{lab} ({desc})")
    return labels_with_desc


"""
Simple utility to check the labels present in a NIfTI segmentation and optionally validate required labels.
Example usage:
    1. Print all labels present in the image
        - ogoCheckLabels segmentation.nii.gz
    2. Check that required labels 2 and 10 are present
        - ogoCheckLabels segmentation.nii.gz --require 2 10
"""
def main():
    parser = argparse.ArgumentParser(
        prog="ogoCheckLabels",
        description="Print labels present in a NIfTI segmentation and optionally validate required labels."
    )
    parser.add_argument("input_image", help="Input image file (*.nii, *.nii.gz)")
    parser.add_argument("--require", type=int, nargs="*", default=[],
                        help="Required labels (e.g. --require 2 10)")
    parser.add_argument("--include-description", action="store_true",
                        help="Include label descriptions in the output")
    args = parser.parse_args()

    ogo.pass_check_if_file_exists(args.input_image)
    ogo.pass_check_file_ending(args.input_image)

    ct = sitk.ReadImage(args.input_image, sitk.sitkUInt8)
    arr = sitk.GetArrayViewFromImage(ct)

    if args.require:
        missing = [lab for lab in args.require if not np.any(arr == lab)]
        if missing:
            if args.include_description:
                missing_with_desc = add_label_descriptions(missing)
                ogo.message("[ERROR] Missing required labels: {}".format(", ".join(missing_with_desc)))
            else:
                ogo.message("[ERROR] Missing required labels: {}".format(missing))
            sys.exit(2)
        ogo.message("[SUCCESS] All required labels are present.")
    else:
        labels = sorted(np.unique(arr).tolist())
        if args.include_description:
            labels_with_desc = add_label_descriptions(labels)
            ogo.message("Labels found: {}".format(", ".join(labels_with_desc)))
        else:
            ogo.message("Labels found: {}".format(labels))

if __name__ == "__main__":
    main()