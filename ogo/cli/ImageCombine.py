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
    validate_integer_range,
    validate_output_file,
    validate_extent,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)

#------------------------------------------------------------------------------+
# Main function
def ImageCombine(input_file1, input_file2, output_file, mode, position, side_position, spacer, filler, set_extent, overwrite, quiet):

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Read input image 1 first to get dimensions if side_position is used
    if not quiet:
        ogo.message(f'Reading input file 1: {input_file1}')
    
    input_image1 = sitk.ReadImage(input_file1)
    
    if not quiet:
        print_image_info(input_file1, input_image1)
    
    # Check that position and side_position are not both set
    if position != [0, 0, 0] and side_position is not None:
        ogo.message('[ERROR] Cannot set both --position and --side_position. Use one or the other.')
        return 1
    
    # If side_position is set, calculate position from image1 dimensions
    if side_position is not None:
        dim = input_image1.GetSize()
        if side_position == 'X':
            position = [dim[0] + spacer, 0, 0]
        elif side_position == 'Y':
            position = [0, dim[1] + spacer, 0]
        elif side_position == 'Z':
            position = [0, 0, dim[2] + spacer]
        if not quiet:
            ogo.message(f'Side position {side_position} with spacer {spacer} set position to: {position}')
    
    # Display parameters
    if not quiet:
        ogo.message(f'Mode: {mode}')
        ogo.message(f'Position: {position}')
        ogo.message(f'Filler value: {filler}')
    
    if not quiet:
        ogo.message(f'Reading input file 2: {input_file2}')
    
    input_image2 = sitk.ReadImage(input_file2)
    
    if not quiet:
        print_image_info(input_file2, input_image2)
    
    # Check that both images have the same spacing
    spacing1 = input_image1.GetSpacing()
    spacing2 = input_image2.GetSpacing()
    
    if spacing1 != spacing2:
        ogo.message('[ERROR] Image spacings do not match.')
        ogo.message(f'       input_file1 spacing: {spacing1}')
        ogo.message(f'       input_file2 spacing: {spacing2}')
        return 1
    
    # Determine pixel type for combined image
    pixel_type1 = input_image1.GetPixelID()
    pixel_type2 = input_image2.GetPixelID()
    
    # Check for unusual pixel types
    if pixel_type1 != sitk.sitkUInt8 and pixel_type1 != sitk.sitkInt16:
        ogo.message(f'[WARNING] input_file1 has unusual pixel type: {sitk.GetPixelIDValueAsString(pixel_type1)}')
    
    if pixel_type2 != sitk.sitkUInt8 and pixel_type2 != sitk.sitkInt16:
        ogo.message(f'[WARNING] input_file2 has unusual pixel type: {sitk.GetPixelIDValueAsString(pixel_type2)}')
    
    # Use UInt8 if both images are UInt8, otherwise use Int16
    if pixel_type1 == sitk.sitkUInt8 and pixel_type2 == sitk.sitkUInt8:
        combined_pixel_type = sitk.sitkUInt8
    else:
        combined_pixel_type = sitk.sitkInt16
    
    if not quiet:
        ogo.message(f'Combined image pixel type: {sitk.GetPixelIDValueAsString(combined_pixel_type)}')
    
    # Validate filler value for UInt8 and Int8 output
    if combined_pixel_type == sitk.sitkUInt8:
        if filler < 0:
            ogo.message(f'[WARNING] Filler value {filler} is out of range for UInt8 (0-255). Setting filler to 0.')
            filler = 0
        elif filler > 255:
            ogo.message(f'[WARNING] Filler value {filler} is out of range for UInt8 (0-255). Setting filler to 255.')
            filler = 255
    elif combined_pixel_type == sitk.sitkInt8:
        if filler < -128:
            ogo.message(f'[WARNING] Filler value {filler} is out of range for Int8 (-128-127). Setting filler to -128.')
            filler = -128
        elif filler > 127:
            ogo.message(f'[WARNING] Filler value {filler} is out of range for Int8 (-128-127). Setting filler to 127.')
            filler = 127
    
    # Get dimensions of both images
    dim1 = input_image1.GetSize()
    dim2 = input_image2.GetSize()
    
    # Calculate extent of combined image
    if set_extent is not None:
        # Use user-provided extent
        minX, maxX, minY, maxY, minZ, maxZ = set_extent
        if not quiet:
            ogo.message(f'Using user-provided extent: X[{minX}:{maxX}], Y[{minY}:{maxY}], Z[{minZ}:{maxZ}]')
        
        # Validate that both images fit within the user-provided extent
        # Image1 is at [0, 0, 0]
        if 0 < minX or dim1[0] > maxX:
            ogo.message(f'[ERROR] Image1 (X extent [0:{dim1[0]}]) does not fit within specified extent X[{minX}:{maxX}]')
            return 1
        if 0 < minY or dim1[1] > maxY:
            ogo.message(f'[ERROR] Image1 (Y extent [0:{dim1[1]}]) does not fit within specified extent Y[{minY}:{maxY}]')
            return 1
        if 0 < minZ or dim1[2] > maxZ:
            ogo.message(f'[ERROR] Image1 (Z extent [0:{dim1[2]}]) does not fit within specified extent Z[{minZ}:{maxZ}]')
            return 1
        
        # Image2 is at position
        if position[0] < minX or (position[0] + dim2[0]) > maxX:
            ogo.message(f'[ERROR] Image2 (X extent [{position[0]}:{position[0] + dim2[0]}]) does not fit within specified extent X[{minX}:{maxX}]')
            return 1
        if position[1] < minY or (position[1] + dim2[1]) > maxY:
            ogo.message(f'[ERROR] Image2 (Y extent [{position[1]}:{position[1] + dim2[1]}]) does not fit within specified extent Y[{minY}:{maxY}]')
            return 1
        if position[2] < minZ or (position[2] + dim2[2]) > maxZ:
            ogo.message(f'[ERROR] Image2 (Z extent [{position[2]}:{position[2] + dim2[2]}]) does not fit within specified extent Z[{minZ}:{maxZ}]')
            return 1
    else:
        # Calculate extent automatically based on image positions
        # Image1 is at [0, 0, 0], Image2 is at position
        minX = min(0, position[0])
        maxX = max(dim1[0], position[0] + dim2[0])
        minY = min(0, position[1])
        maxY = max(dim1[1], position[1] + dim2[1])
        minZ = min(0, position[2])
        maxZ = max(dim1[2], position[2] + dim2[2])
    
    # Dimensions of combined image
    combined_dim = [maxX - minX, maxY - minY, maxZ - minZ]
    
    if not quiet:
        ogo.message(f'Combined image dimensions: {combined_dim}')
        ogo.message(f'Combined image extent (voxels): X[{minX}:{maxX}], Y[{minY}:{maxY}], Z[{minZ}:{maxZ}]')
        
        # Calculate physical extent using image1's coordinate system
        phys_min = input_image1.TransformContinuousIndexToPhysicalPoint([float(minX), float(minY), float(minZ)])
        phys_max = input_image1.TransformContinuousIndexToPhysicalPoint([float(maxX), float(maxY), float(maxZ)])
        ogo.message(f'Combined image extent (physical): X[{phys_min[0]:.2f}:{phys_max[0]:.2f}], Y[{phys_min[1]:.2f}:{phys_max[1]:.2f}], Z[{phys_min[2]:.2f}:{phys_max[2]:.2f}] mm')
    
    # Create combined image initialized with filler value
    combined_image = sitk.Image(combined_dim, combined_pixel_type)
    combined_image.SetSpacing(input_image1.GetSpacing())
    
    # Adjust origin to account for extent that may extend beyond image1
    index_offset = [float(minX), float(minY), float(minZ)]
    new_origin = input_image1.TransformContinuousIndexToPhysicalPoint(index_offset)
    combined_image.SetOrigin(new_origin)
    
    combined_image.SetDirection(input_image1.GetDirection())
    
    # Initialize with filler value
    combined_image = combined_image + filler
    
    if not quiet:
        ogo.message(f'Created combined image initialized with filler value: {filler}')
    
    # Paste image1 into combined_image
    # Image1 is at [0, 0, 0] in original coordinate system
    # In combined image, it's at position [0-minX, 0-minY, 0-minZ]
    paste_position1 = [0 - minX, 0 - minY, 0 - minZ]
    
    if not quiet:
        ogo.message(f'Pasting image1 at position: {paste_position1}')
    
    combined_image = sitk.Paste(
        combined_image,
        input_image1,
        input_image1.GetSize(),
        [0, 0, 0],
        paste_position1
    )
    
    # Paste image2 into combined_image based on mode
    # Image2 is at position in original coordinate system
    # In combined image, it's at position [position[0]-minX, position[1]-minY, position[2]-minZ]
    paste_position2 = [position[0] - minX, position[1] - minY, position[2] - minZ]
    
    if not quiet:
        ogo.message(f'Pasting image2 at position: {paste_position2}')
    
    if mode == 'overlay':
        # Overlay mode: Paste image2, overwriting existing values except where image2 is 0
        # Extract the region where image2 will go
        region = sitk.RegionOfInterest(
            combined_image,
            input_image2.GetSize(),
            paste_position2
        )
        # Cast image2 to match combined_image pixel type
        image2_casted = sitk.Cast(input_image2, combined_pixel_type)
        # Copy physical space information from region to match
        image2_casted.CopyInformation(region)
        
        # Create mask where image2 is non-zero
        mask = sitk.Cast(image2_casted != 0, combined_pixel_type)
        
        # Combine: use image2 where mask is non-zero, otherwise keep region
        result = mask * image2_casted + (1 - mask) * region
        
        # Paste the result back
        combined_image = sitk.Paste(
            combined_image,
            result,
            result.GetSize(),
            [0, 0, 0],
            paste_position2
        )
    elif mode == 'add':
        # Add mode: Add image2 values to existing combined_image values
        # Extract the region where image2 will go
        region = sitk.RegionOfInterest(
            combined_image,
            input_image2.GetSize(),
            paste_position2
        )
        # Cast image2 to match combined_image pixel type for arithmetic operation
        image2_casted = sitk.Cast(input_image2, combined_pixel_type)
        # Copy physical space information from region to match for arithmetic operation
        image2_casted.CopyInformation(region)
        # Add image2 to this region
        region = sitk.Add(region, image2_casted)
        # Paste the result back
        combined_image = sitk.Paste(
            combined_image,
            region,
            region.GetSize(),
            [0, 0, 0],
            paste_position2
        )
    elif mode == 'subtract':
        # Subtract mode: Subtract image2 values from existing combined_image values
        # Extract the region where image2 will go
        region = sitk.RegionOfInterest(
            combined_image,
            input_image2.GetSize(),
            paste_position2
        )
        # Cast image2 to match combined_image pixel type for arithmetic operation
        image2_casted = sitk.Cast(input_image2, combined_pixel_type)
        # Copy physical space information from region to match for arithmetic operation
        image2_casted.CopyInformation(region)
        # Subtract image2 from this region
        region = sitk.Subtract(region, image2_casted)
        # Paste the result back
        combined_image = sitk.Paste(
            combined_image,
            region,
            region.GetSize(),
            [0, 0, 0],
            paste_position2
        )
    
    if not quiet:
        ogo.message(f'Writing output file: {output_file}')
    
    sitk.WriteImage(combined_image, output_file)
    
    if not quiet:
        print_image_info(output_file, combined_image)
        
        # Display final extent of output image in both voxels and physical coordinates
        ogo.message(f'Output image extent (voxels): X[{minX}:{maxX}], Y[{minY}:{maxY}], Z[{minZ}:{maxZ}]')
        
        output_origin = combined_image.GetOrigin()
        output_size = combined_image.GetSize()
        output_spacing = combined_image.GetSpacing()
        output_phys_max = combined_image.TransformContinuousIndexToPhysicalPoint(
            [float(output_size[0]), float(output_size[1]), float(output_size[2])]
        )
        ogo.message(f'Output image extent (physical): X[{output_origin[0]:.2f}:{output_phys_max[0]:.2f}], Y[{output_origin[1]:.2f}:{output_phys_max[1]:.2f}], Z[{output_origin[2]:.2f}:{output_phys_max[2]:.2f}] mm')
    
    if not quiet:
        ogo.message('[DONE]')
    
    return 0

def main():
    # Setup description
    description='''
Utility to combine two NIFTI image files with flexible positioning and blending modes.
This performs the combination in voxel units, not physical space (i.e., Origin).

This tool takes two input NIFTI images and combines them according to the specified mode.
It is not required that the two input images have compatible dimensions. The output image 
will have dimensions that cover the extent of the two input images. The filler variable 
will be used in regions where neither image provides data within the extent.

The mode determines how overlapping voxels are handled:
  - add: Add voxel values where images overlap (image 1 + image 2)
  - subtract: Subtract voxel values where images overlap (image 1 - image 2)
  - overlay: Prioritize image 2 voxels (overwrite image 1 voxels)

The position parameter defines the relative position of input_file2 with respect to 
input_file1 in voxel coordinates. Alternatively, use --side_position for convenient 
positioning along dimension boundaries of input_file1 (X, Y, or Z axis). The --spacer
parameter can be used with --side_position to add additional spacing (in voxels) between
the images.

Optionally, use --set_extent to manually specify the extent of the output combined image.
'''
    epilog='''
Example calls: 
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz --mode overlay --position 10 20 30
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz -m add --filler 255
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz -m subtract -p -5 0 10 -ow -q
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz --side_position X
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz -s Z -m overlay -ow
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz -s X --spacer 10
ogoImageCombine input1.nii.gz input2.nii.gz output.nii.gz -e -10 100 0 200 0 150 -m add
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageCombine",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file1', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='First input image file (*.nii, *.nii.gz)')
    parser.add_argument('input_file2', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Second input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.nii', '.nii.gz']), 
                        help='Output combined image file (*.nii, *.nii.gz)')
    parser.add_argument('--mode', '-m', default='overlay', choices=['add', 'subtract', 'overlay'],
                        help='Mode for combining overlapping voxels (default: %(default)s)')
    parser.add_argument('--position', '-p', type=int, nargs=3, default=[0, 0, 0], metavar='INT',
                        help='Relative position of input_file2 [x, y, z] in voxels (default: 0 0 0)')
    parser.add_argument('--side_position', '-s', default=None, choices=['X', 'Y', 'Z'],
                        help='Position input_file2 at dimension boundary (X, Y, or Z) of input_file1. Alternative to --position')
    parser.add_argument('--spacer', '-sp', type=validate_integer_min(0), default=0, metavar='INT',
                        help='Additional spacing (in voxels) when using --side_position (must be >= 0, default: %(default)s)')
    parser.add_argument('--filler', '-f', type=validate_integer_range(-32768, 32767), default=0, metavar='VALUE',
                        help='Filler value for regions without data (default: %(default)s)')
    parser.add_argument('--set_extent', '-e', type=int, nargs=6, default=None, metavar=('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'),
                        help='Manually set extent of combined image [xmin xmax ymin ymax zmin zmax] in voxels')
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
    if args.input_file1 is None or args.input_file2 is None or args.output_file is None:
        parser.error('input_file1, input_file2, and output_file are required')
    
    # Validate set_extent if provided
    if args.set_extent is not None:
        try:
            args.set_extent = validate_extent(args.set_extent)
        except argparse.ArgumentTypeError as e:
            parser.error(str(e))
    
    if not args.quiet:
        ogo.message(echo_arguments('ImageCombine', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = ImageCombine(**run_args)
    return result

if __name__ == '__main__':
    main()
