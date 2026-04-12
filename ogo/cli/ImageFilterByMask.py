# /------------------------------------------------------------------------------+
# | 11-APR-2026                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import SimpleITK as sitk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
from ogo.util.cli_tools import (
    check_overwrite,
    validate_integer_min,
    validate_input_file,
    validate_output_file,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)

#------------------------------------------------------------------------------+
# Main function
def ImageFilterByMask(input_file, mask_file, output_file, mask_dilate, mask_erode, morphological, overwrite, quiet):

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Read input image
    if not quiet:
        ogo.message(f'Reading input file: {input_file}')
    
    input_image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, input_image)
    
    # Check pixel type
    if input_image.GetPixelID() != sitk.sitkInt16:
        ogo.message(f'[ERROR] Input file must be Int16, got {sitk.GetPixelIDValueAsString(input_image.GetPixelID())}')
        return 1
    
    # Read mask
    if not quiet:
        ogo.message(f'Reading mask file: {mask_file}')
    
    mask_image = sitk.ReadImage(mask_file)
    
    if not quiet:
        print_image_info(mask_file, mask_image)
    
    # Check mask pixel type
    if mask_image.GetPixelID() != sitk.sitkUInt8:
        ogo.message(f'[ERROR] Mask file must be UInt8, got {sitk.GetPixelIDValueAsString(mask_image.GetPixelID())}')
        return 1
    
    # Check that images have matching dimensions
    if input_image.GetSize() != mask_image.GetSize():
        ogo.message('[ERROR] Input and mask dimensions do not match.')
        ogo.message(f'       input_file: {input_image.GetSize()}')
        ogo.message(f'       mask_file: {mask_image.GetSize()}')
        return 1
    
    # Convert mask to binary (0 or 1) for morphological operations
    processed_mask = sitk.Cast(mask_image != 0, sitk.sitkUInt8)
    
    # Warn if both dilate and erode are specified without explicit morphological parameter
    if mask_dilate > 0 and mask_erode > 0 and morphological is None:
        ogo.message('[WARNING] Both dilate and erode specified. Dilation will be applied first, then erosion.')
        ogo.message('          This is equivalent to a morphological closing operation.')
        ogo.message('          To explicitly control operation order, use --morphological opening or --morphological closing.')
    
    # Warn if morphological operation selected but no dilate/erode specified
    if morphological in ['opening', 'closing'] and mask_dilate == 0 and mask_erode == 0:
        ogo.message(f'[WARNING] Morphological {morphological} selected but neither --mask_dilate nor --mask_erode specified.')
        ogo.message('          No morphological operations will be applied to the mask.')
    
    # Determine operation order based on morphological parameter
    if morphological == 'opening':
        # Opening: erode then dilate (removes small objects, smooths boundaries from outside)
        if not quiet:
            ogo.message(f'Using morphological opening (erode: {mask_erode} voxels, dilate: {mask_dilate} voxels)')
        
        if mask_erode > 0:
            if not quiet:
                ogo.message(f'Eroding mask by {mask_erode} voxels (opening: step 1/2)')
            erode_filter = sitk.BinaryErodeImageFilter()
            erode_filter.SetKernelRadius(mask_erode)
            processed_mask = erode_filter.Execute(processed_mask)
        
        if mask_dilate > 0:
            if not quiet:
                ogo.message(f'Dilating mask by {mask_dilate} voxels (opening: step 2/2)')
            dilate_filter = sitk.BinaryDilateImageFilter()
            dilate_filter.SetKernelRadius(mask_dilate)
            processed_mask = dilate_filter.Execute(processed_mask)
    
    elif morphological == 'closing':
        # Closing: dilate then erode (fills small holes, smooths boundaries from inside)
        if not quiet:
            ogo.message(f'Using morphological closing (dilate: {mask_dilate} voxels, erode: {mask_erode} voxels)')
        
        if mask_dilate > 0:
            if not quiet:
                ogo.message(f'Dilating mask by {mask_dilate} voxels (closing: step 1/2)')
            dilate_filter = sitk.BinaryDilateImageFilter()
            dilate_filter.SetKernelRadius(mask_dilate)
            processed_mask = dilate_filter.Execute(processed_mask)
        
        if mask_erode > 0:
            if not quiet:
                ogo.message(f'Eroding mask by {mask_erode} voxels (closing: step 2/2)')
            erode_filter = sitk.BinaryErodeImageFilter()
            erode_filter.SetKernelRadius(mask_erode)
            processed_mask = erode_filter.Execute(processed_mask)
    
    else:
        # Default behavior: dilate then erode
        if mask_dilate > 0:
            if not quiet:
                ogo.message(f'Dilating mask by {mask_dilate} voxels')
            dilate_filter = sitk.BinaryDilateImageFilter()
            dilate_filter.SetKernelRadius(mask_dilate)
            processed_mask = dilate_filter.Execute(processed_mask)
        
        if mask_erode > 0:
            if not quiet:
                ogo.message(f'Eroding mask by {mask_erode} voxels')
            erode_filter = sitk.BinaryErodeImageFilter()
            erode_filter.SetKernelRadius(mask_erode)
            processed_mask = erode_filter.Execute(processed_mask)
    
    # Create binary mask as Int16: 1 where mask is non-zero, 0 where mask is zero
    binary_mask = sitk.Cast(processed_mask, sitk.sitkInt16)
    
    # Apply mask: where mask=0, output=0; where mask!=0, output=input value
    if not quiet:
        ogo.message('Applying mask to input image')
    
    output_image = input_image * binary_mask
    
    # Ensure output is Int16
    output_image = sitk.Cast(output_image, sitk.sitkInt16)
    
    # Copy metadata from input image
    output_image.CopyInformation(input_image)
    
    # Write output
    if not quiet:
        ogo.message(f'Writing output file: {output_file}')
    
    sitk.WriteImage(output_image, output_file)
    
    if not quiet:
        print_image_info(output_file, output_image)
        ogo.message('[DONE]')
    
    return 0

def main():
    # Setup description
    description='''
Utility to apply a binary mask to an input image.

This tool applies a mask to set voxels to zero where the mask is zero, and passes 
through the input image values where the mask is non-zero. The input image must be 
Int16 and the mask must be UInt8. The output will always be Int16.

Optional morphological operations (dilation and erosion) can be applied to the 
mask before applying it to the input image.

Morphological operation types:
  - opening: Erode then dilate (removes small objects, smooths from outside)
  - closing: Dilate then erode (fills small holes, smooths from inside)
  - Default (no --morphological): Dilate then erode if both specified
'''
    epilog='''
Example calls: 
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz --mask_dilate 2
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz --mask_erode 1 -ow
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz -d 3 -e 1 -q
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz -d 2 -e 2 --morphological opening
ogoImageMask input.nii.gz mask.nii.gz output.nii.gz -d 3 -e 3 -mo closing
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageMask",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Input image file (*.nii, *.nii.gz) - must be Int16')
    parser.add_argument('mask_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Mask image file (*.nii, *.nii.gz) - must be UInt8')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.nii', '.nii.gz']), 
                        help='Output masked image file (*.nii, *.nii.gz)')
    parser.add_argument('--mask_dilate', '-d', type=validate_integer_min(0), default=0, metavar='VOXELS',
                        help='Dilate mask by specified number of voxels (default: %(default)s)')
    parser.add_argument('--mask_erode', '-e', type=validate_integer_min(0), default=0, metavar='VOXELS',
                        help='Erode mask by specified number of voxels (default: %(default)s)')
    parser.add_argument('--morphological', '-mo', choices=['opening', 'closing'], default=None, metavar='TYPE',
                        help='Morphological operation type: opening (erode then dilate) or closing (dilate then erode)')
    parser.add_argument('--overwrite','-ow', action='store_true', help='Overwrite output without asking')
    parser.add_argument('--quiet','-q', action='store_true', help='Suppress informational output to stdout')

    print()

    # Parse and display
    args = parser.parse_args()
    
    # Handle --json flag
    if args.json is not None:
        if args.json == '':
            # --json with no argument: output to stdout
            save_config_to_json(args, None)
            return 0
        else:
            # --json file.json: load from file
            if not args.json.endswith('.json'):
                parser.error('JSON file must have .json extension')
            json_config = load_config_from_json(args.json)
            # Update args with JSON values (command line args override JSON)
            for key, value in json_config.items():
                if hasattr(args, key):
                    # Only update if argument was not explicitly provided on command line
                    if getattr(args, key) == parser.get_default(key):
                        setattr(args, key, value)
    
    # Check that required positional arguments are provided
    if args.input_file is None or args.mask_file is None or args.output_file is None:
        parser.error('input_file, mask_file, and output_file are required')
    
    if not args.quiet:
        ogo.message(echo_arguments('ImageFilterByMask', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = ImageFilterByMask(**run_args)
    return result

if __name__ == '__main__':
    main()
