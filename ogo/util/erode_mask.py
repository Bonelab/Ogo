import SimpleITK as sitk
import argparse
import sys
from ogo.util.echo_arguments import echo_arguments

#function that will erode a vertebral body mask by 2mm in each dimension (x,y,z) except if the slice thickness is too big
#if the slice thickness is >= 2mm, then we will erode by 1 voxel instead of 2mm to avoid losing too much of the mask
def erode_mask(mask_path, output_path, label):
    
    # Read the mask
    mask = sitk.ReadImage(mask_path)
    spacing = mask.GetSpacing()
    
    print(f"Mask shape: {mask.GetSize()}")
    print(f"Spacing (voxel size): {spacing}")
    
    # determine erosion radius in voxels for each dimension separately
    erosion_radius = []
    for idx, s in enumerate(spacing):
        if s >= 2.0:
            erosion_radius.append(1)
        else:
            erosion_radius.append(max(1, round(2.0 / s)))
    
    print(f"Erosion radius in voxels per dimension: {erosion_radius}")
    
    binary_mask = sitk.BinaryThreshold(mask, lowerThreshold=label, upperThreshold=label)
    

    print(f"Erosion radius: {erosion_radius}")

    eroded = sitk.BinaryErode(binary_mask, erosion_radius)

    # make sure its the original mask label (20 for L1)
    eroded = sitk.Cast(eroded, mask.GetPixelID())
    eroded = sitk.Multiply(eroded, label)
    eroded.CopyInformation(mask)

    sitk.WriteImage(eroded, output_path)
    
    print(f"Eroded mask saved to: {output_path}")

    return eroded


def main():
    parser = argparse.ArgumentParser(
        description="Erode a NIFTI mask by 2mm (or 1 voxel for thick slices)"
    )
    parser.add_argument("mask_path", help="Path to input mask NIFTI file")
    parser.add_argument("output_path", help="Path to save eroded mask")
    parser.add_argument("--label", type=int, default=20, help="Label value to erode (default: 20)")
    
    args = parser.parse_args()
    
    try:
        erode_mask(
            mask_path=args.mask_path,
            output_path=args.output_path,
            label=args.label
        )
        print(f"Success: Eroded mask saved to {args.output_path}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
