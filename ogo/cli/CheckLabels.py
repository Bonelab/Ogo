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

def get_labels(ct):
    filt = sitk.LabelShapeStatisticsImageFilter()
    filt.Execute(ct)
    return list(filt.GetLabels())

def main():
    parser = argparse.ArgumentParser(
        prog="ogoCheckLabels",
        description="Print labels present in a NIfTI segmentation and optionally validate required labels."
    )
    parser.add_argument("input_image", help="Input image file (*.nii, *.nii.gz)")
    parser.add_argument("--require", type=int, nargs="*", default=[],
                        help="Required labels (e.g. --require 2 10)")
    args = parser.parse_args()

    if not os.path.isfile(args.input_image):
        sys.exit(f'[ERROR] Cannot find file "{args.input_image}"')
    if not (args.input_image.lower().endswith(".nii") or args.input_image.lower().endswith(".nii.gz")):
        sys.exit(f'[ERROR] Input must be type NIFTI file: "{args.input_image}"')

    ct = sitk.ReadImage(args.input_image, sitk.sitkUInt8)
    labels = sorted(get_labels(ct))
    print("Labels found:", labels)

    if args.require:
        missing = [lab for lab in args.require if lab not in labels]
        if missing:
            print("[ERROR] Missing required labels:", missing)
            sys.exit(2)
        print("[SUCCESS] All required labels are present.")

if __name__ == "__main__":
    main()