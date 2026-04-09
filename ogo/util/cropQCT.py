import SimpleITK as sitk
import numpy as np
import argparse

from ogo.util.echo_arguments import echo_arguments

def load_image_and_mask(image_path, mask_path):
    image = sitk.ReadImage(image_path)
    mask = sitk.ReadImage(mask_path)
    return image, mask

def find_mask_bounds(mask_sitk):
    mask_array = sitk.GetArrayFromImage(mask_sitk)
    non_zero_indices = mask_array.nonzero()
    min_x, max_x = non_zero_indices[0].min(), non_zero_indices[0].max()
    min_y, max_y = non_zero_indices[1].min(), non_zero_indices[1].max()
    min_z, max_z = non_zero_indices[2].min(), non_zero_indices[2].max()
    return (min_x, max_x), (min_y, max_y), (min_z, max_z)

def crop_image(image_sitk, mask_sitk, buffer=30):

    #Crops an image and its corresponding mask image with a buffer on the outside so that it is not exactly down to the mask.
    
    (min_x, max_x), (min_y, max_y), (min_z, max_z) = find_mask_bounds(mask_sitk)
    
    min_x = max(min_x - buffer, 0)
    max_x = min(max_x + buffer, mask_sitk.GetSize()[2] - 1)
    min_y = max(min_y - buffer, 0)
    max_y = min(max_y + buffer, mask_sitk.GetSize()[1] - 1)
    min_z = max(min_z - buffer, 0)
    max_z = min(max_z + buffer, mask_sitk.GetSize()[0] - 1)
    
    extract_size = [int(max_z - min_z + 1), int(max_y - min_y + 1), int(max_x - min_x + 1)]
    extract_index = [int(min_z), int(min_y), int(min_x)]
    
    extractor = sitk.RegionOfInterestImageFilter()
    extractor.SetSize(extract_size)
    extractor.SetIndex(extract_index)
    
    cropped_image = extractor.Execute(image_sitk)
    cropped_mask = extractor.Execute(mask_sitk)
    
    return cropped_image, cropped_mask


def main():
    description = '''
    This is a helper script for reading in an image and a mask and outputting a cropped image + mask based on specified bounds.
    This is mainly for creating cropped QCT images that can be used with the extracted DXA femoral neck and total hip ROIs.



'''
    epilog = '''
    Example calls:

    python cropQCT.py image_k2hpo4.nii mask.nii.gz cropped_image.nii.gz cropped_mask.nii.gz --buffer 30

    '''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="cropQCT",
        description=description,
        epilog=epilog
    )

    parser.add_argument('image_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('mask_filename', help='Input image mask file (*.nii, *.nii.gz)')
    parser.add_argument('cropped_image',
                        help='Path to where the output cropped image will be output')
    parser.add_argument('cropped_mask', help='Output file name for cropped mask')
    parser.add_argument('--buffer', type=int, default=30, help='Buffer size in pixels (default: %(default)s)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    print()

    # Parse and display
    args = parser.parse_args()

    if args.cropped_image:
        print(echo_arguments('cropQCT', vars(args)))

    # Run program
    image, mask = load_image_and_mask(args.image_filename, args.mask_filename)
    cropped_image, cropped_mask = crop_image(image, mask, buffer=args.buffer)
    sitk.WriteImage(cropped_image, args.cropped_image)


if __name__ == '__main__':
    main()