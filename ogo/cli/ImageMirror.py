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
    validate_input_file,
    validate_output_file,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)

#------------------------------------------------------------------------------+
# Main function
def ImageMirror(input_file, output_file, mirror_plane, flip_about_origin, overwrite, quiet):

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Read input image
    if not quiet:
        ogo.message(f'Reading input file: {input_file}')
    
    input_image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, input_image)
    
    # Determine which axis to flip based on mirror plane
    flip_axes = [False, False, False]
    
    if mirror_plane == 'X':
        flip_axes[0] = True  # Flip X axis
        axis_name = 'X'
    elif mirror_plane == 'Y':
        flip_axes[1] = True  # Flip Y axis
        axis_name = 'Y'
    elif mirror_plane == 'Z':
        flip_axes[2] = True  # Flip Z axis
        axis_name = 'Z'
    
    if not quiet:
        ogo.message(f'Flipping image along {axis_name} axis')
        if flip_about_origin:
            ogo.message(f'Flip about origin: True (flipping across coordinate system origin)')
        else:
            ogo.message(f'Flip about origin: False (flipping within coordinate grid)')
    
    # Apply flip filter
    # Note: SimpleITK automatically updates origin and direction matrix
    # to preserve physical space. Voxel indices are reversed.
    output_image = sitk.Flip(input_image, flip_axes, flip_about_origin)
    
    if not quiet:
        ogo.message(f'Voxel order reversed along {axis_name} axis')
        if not flip_about_origin:
            ogo.message(f'Physical location of voxels preserved')
    
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
Mirror an image by flipping along a specified axis.

This tool flips an image along one of the three principal axes:
  - X: Flip along the X axis
  - Y: Flip along the Y axis
  - Z: Flip along the Z axis

The mirroring operation reverses the order of voxels along the specified axis.

IMPORTANT: SimpleITK automatically updates the origin and direction matrix of the 
image so that the physical location of the pixels remains consistent. If you view 
the flipped image in a physical-space-aware tool like ITK-SNAP or 3D Slicer, it 
may look the same as the original, but the voxel indices (i, j, k) will be reversed.

Flip modes:
  - flip_about_origin=False (default): Flips within the coordinate grid,
    preserving physical location of voxels
  - flip_about_origin=True: Flips across the coordinate system's origin,
    which may change the physical location
'''
    epilog='''
Example calls: 
ogoImageMirror input.nii.gz output.nii.gz X
ogoImageMirror input.nii.gz output.nii.gz Y -ow
ogoImageMirror input.nii.gz output.nii.gz Z --quiet
ogoImageMirror input.nii.gz mirrored.nii.gz X --flip_about_origin
ogoImageMirror input.nii.gz mirrored.nii.gz Y -fao -q
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageMirror",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.nii', '.nii.gz']), 
                        help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('mirror_plane', nargs='?', choices=['X', 'Y', 'Z'],
                        help='Axis to flip: X, Y, or Z')
    parser.add_argument('--flip_about_origin', '-fao', action='store_true', default=False,
                        help='Flip across coordinate system origin (default: flip within coordinate grid)')
    parser.add_argument('--overwrite', '-ow', action='store_true', help='Overwrite output without asking')
    parser.add_argument('--quiet', '-q', action='store_true', help='Suppress informational output to stdout')

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
    if args.input_file is None or args.output_file is None or args.mirror_plane is None:
        parser.error('input_file, output_file, and mirror_plane are required')
    
    if not args.quiet:
        ogo.message(echo_arguments('ImageMirror', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = ImageMirror(**run_args)
    return result

if __name__ == '__main__':
    main()
