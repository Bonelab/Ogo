import os
import argparse
import SimpleITK as sitk

def main():
    """
    Converts one or more NIfTI files to RAS orientation using SimpleITK.
    Accepts either a single file or a folder of NIfTI files.
    """
    parser = argparse.ArgumentParser(description="Convert NIfTI file(s) to RAS orientation using SimpleITK.")
    parser.add_argument("input_path", type=str, help="Path to a NIfTI file or folder containing NIfTI files.")
    parser.add_argument("-o", "--output_folder", type=str, default=None, help="Optional path to save RAS images.")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress verbose output.")

    args = parser.parse_args()
    input_path = args.input_path
    output_folder = args.output_folder
    verbose = not args.quiet

    def get_orientation(img):
        direction = img.GetDirection()
        return sitk.DICOMOrientImageFilter().GetOrientationFromDirectionCosines(direction)

    orienter = sitk.DICOMOrientImageFilter()
    orienter.SetDesiredCoordinateOrientation("RAS")

    if os.path.isfile(input_path) and (input_path.endswith(".nii") or input_path.endswith(".nii.gz")):
        # Single file case
        image = sitk.ReadImage(input_path)
        original_orientation = get_orientation(image)
        image_ras = orienter.Execute(image)
        final_orientation = get_orientation(image_ras)

        if verbose:
            print(f"{os.path.basename(input_path)}")
            print(f"  Original orientation: {original_orientation}")
            print(f"  Final orientation:    {final_orientation}")

        output_dir = output_folder or os.path.dirname(input_path)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, os.path.basename(input_path))
        sitk.WriteImage(image_ras, output_path)

    elif os.path.isdir(input_path):
        # Folder case
        output_dir = output_folder or input_path
        os.makedirs(output_dir, exist_ok=True)

        for fname in os.listdir(input_path):
            if fname.endswith(".nii") or fname.endswith(".nii.gz"):
                file_path = os.path.join(input_path, fname)
                image = sitk.ReadImage(file_path)

                original_orientation = get_orientation(image)
                image_ras = orienter.Execute(image)
                final_orientation = get_orientation(image_ras)

                if verbose:
                    print(f"{fname}")
                    print(f"  Original orientation: {original_orientation}")
                    print(f"  Final orientation:    {final_orientation}")

                output_path = os.path.join(output_dir, fname)
                sitk.WriteImage(image_ras, output_path)
    else:
        raise ValueError("Input path must be a NIfTI file (.nii or .nii.gz) or a folder containing such files.")

if __name__ == "__main__":
    main()
