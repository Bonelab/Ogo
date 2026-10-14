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
def ImageChangeExtent(input_file, output_file, offset, extent, filler, overwrite, quiet):

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Check that offset and extent are mutually exclusive
    if offset is not None and extent is not None:
        ogo.message('[ERROR] Cannot specify both --offset and --extent.')
        ogo.message('        Please specify one or the other.')
        return 1
    
    # Check that at least one is specified
    if offset is None and extent is None:
        ogo.message('[ERROR] Must specify either --offset or --extent.')
        return 1
    
    # Read input image
    if not quiet:
        ogo.message(f'Reading input file: {input_file}')
    
    input_image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, input_image)
    
    # Get image properties
    input_size = input_image.GetSize()
    input_origin = input_image.GetOrigin()
    input_spacing = input_image.GetSpacing()
    pixel_type = input_image.GetPixelID()
    
    # Process based on offset or extent
    if offset is not None:
        # Validate offset dimensions
        if len(offset) != 3:
            ogo.message(f'[ERROR] Offset must have exactly 3 values, got {len(offset)}')
            return 1
        
        ox, oy, oz = offset
        sx, sy, sz = input_size
        
        # Validate negative offsets don't equal or exceed half dimension
        if ox < 0 and abs(ox) >= sx // 2:
            ogo.message(f'[ERROR] Negative X offset ({ox}) equals or exceeds half the X dimension ({sx // 2})')
            return 1
        if oy < 0 and abs(oy) >= sy // 2:
            ogo.message(f'[ERROR] Negative Y offset ({oy}) equals or exceeds half the Y dimension ({sy // 2})')
            return 1
        if oz < 0 and abs(oz) >= sz // 2:
            ogo.message(f'[ERROR] Negative Z offset ({oz}) equals or exceeds half the Z dimension ({sz // 2})')
            return 1
        
        # Calculate new extent based on offset
        # Offset is applied symmetrically on both sides of each dimension
        xmin, ymin, zmin = -ox, -oy, -oz
        xmax = sx + ox
        ymax = sy + oy
        zmax = sz + oz
        
        if not quiet:
            if ox >= 0 and oy >= 0 and oz >= 0:
                ogo.message(f'Expanding image by offset: [{ox}, {oy}, {oz}]')
            elif ox <= 0 and oy <= 0 and oz <= 0:
                ogo.message(f'Cropping image by offset: [{ox}, {oy}, {oz}]')
            else:
                ogo.message(f'Adjusting image by offset: [{ox}, {oy}, {oz}]')
            ogo.message(f'New extent (voxels): [{xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}]')
    
    else:  # extent is not None
        # Validate extent
        if len(extent) != 6:
            ogo.message(f'[ERROR] Extent must have exactly 6 values, got {len(extent)}')
            return 1
        
        xmin, xmax, ymin, ymax, zmin, zmax = extent
        
        # Validate min < max for each dimension
        if xmin >= xmax:
            ogo.message(f'[ERROR] xmin ({xmin}) must be less than xmax ({xmax})')
            return 1
        if ymin >= ymax:
            ogo.message(f'[ERROR] ymin ({ymin}) must be less than ymax ({ymax})')
            return 1
        if zmin >= zmax:
            ogo.message(f'[ERROR] zmin ({zmin}) must be less than zmax ({zmax})')
            return 1
        
        if not quiet:
            ogo.message(f'Using explicit extent: [{xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}]')
    
    # Calculate new size
    new_size = [xmax - xmin, ymax - ymin, zmax - zmin]
    
    # Create output image filled with filler value
    if not quiet:
        ogo.message(f'Creating output image with size: {new_size}')
        ogo.message(f'Filler value: {filler}')
    
    output_image = sitk.Image(new_size, pixel_type)
    output_image.SetSpacing(input_spacing)
    output_image.SetDirection(input_image.GetDirection())
    
    # Calculate new origin
    # Origin shifts when we change the extent
    index_offset = [float(xmin), float(ymin), float(zmin)]
    new_origin = input_image.TransformContinuousIndexToPhysicalPoint(index_offset)
    output_image.SetOrigin(new_origin)
    
    # Fill with filler value
    if filler != 0:
        output_image = output_image + filler
    
    # Determine the region to paste from input image
    # Calculate overlap between input image extent [0, sx) x [0, sy) x [0, sz)
    # and new extent [xmin, xmax) x [ymin, ymax) x [zmin, zmax)
    
    paste_dest_start = [max(0, -xmin), max(0, -ymin), max(0, -zmin)]
    crop_source_start = [max(0, xmin), max(0, ymin), max(0, zmin)]
    
    overlap_size = [
        min(input_size[0], xmax) - max(0, xmin),
        min(input_size[1], ymax) - max(0, ymin),
        min(input_size[2], zmax) - max(0, zmin)
    ]
    
    # Check if there's any overlap
    if overlap_size[0] <= 0 or overlap_size[1] <= 0 or overlap_size[2] <= 0:
        ogo.message('[WARNING] No overlap between input image and new extent.')
        ogo.message('          Output will contain only filler values.')
    else:
        # Extract the region from input image
        if not quiet:
            ogo.message(f'Extracting region from input: start={crop_source_start}, size={overlap_size}')
            ogo.message(f'Pasting into output at: {paste_dest_start}')
        
        region = sitk.RegionOfInterest(input_image, overlap_size, crop_source_start)
        
        # Paste into output image
        output_image = sitk.Paste(output_image, region, region.GetSize(), [0, 0, 0], paste_dest_start)
    
    # Display extent information
    if not quiet:
        # Input image extent
        input_extent_voxels = [0, input_size[0], 0, input_size[1], 0, input_size[2]]
        input_phys_min = input_image.TransformContinuousIndexToPhysicalPoint([0.0, 0.0, 0.0])
        input_phys_max = input_image.TransformContinuousIndexToPhysicalPoint([float(input_size[0]), float(input_size[1]), float(input_size[2])])
        
        ogo.message(f'Input image extent:')
        ogo.message(f'  Voxels:   [{input_extent_voxels[0]}, {input_extent_voxels[1]}, {input_extent_voxels[2]}, {input_extent_voxels[3]}, {input_extent_voxels[4]}, {input_extent_voxels[5]}]')
        ogo.message(f'  Physical: [{input_phys_min[0]:.2f}, {input_phys_max[0]:.2f}, {input_phys_min[1]:.2f}, {input_phys_max[1]:.2f}, {input_phys_min[2]:.2f}, {input_phys_max[2]:.2f}] mm')
        
        # Output image extent - use input image transformation with extent coordinates
        output_phys_min = input_image.TransformContinuousIndexToPhysicalPoint([float(xmin), float(ymin), float(zmin)])
        output_phys_max = input_image.TransformContinuousIndexToPhysicalPoint([float(xmax), float(ymax), float(zmax)])
        
        ogo.message(f'Output image extent:')
        ogo.message(f'  Voxels:   [{xmin}, {xmax}, {ymin}, {ymax}, {zmin}, {zmax}]')
        ogo.message(f'  Physical: [{output_phys_min[0]:.2f}, {output_phys_max[0]:.2f}, {output_phys_min[1]:.2f}, {output_phys_max[1]:.2f}, {output_phys_min[2]:.2f}, {output_phys_max[2]:.2f}] mm')
        
        # Calculate voxel statistics
        total_input_voxels = input_size[0] * input_size[1] * input_size[2]
        total_output_voxels = new_size[0] * new_size[1] * new_size[2]
        
        if overlap_size[0] > 0 and overlap_size[1] > 0 and overlap_size[2] > 0:
            pasted_voxels = overlap_size[0] * overlap_size[1] * overlap_size[2]
            percent_pasted = (pasted_voxels / total_input_voxels) * 100.0
        else:
            pasted_voxels = 0
            percent_pasted = 0.0
        
        ogo.message(f'Voxel statistics:')
        ogo.message(f'  Input voxels:  {total_input_voxels:,}')
        ogo.message(f'  Output voxels: {total_output_voxels:,}')
        ogo.message(f'  Pasted voxels: {pasted_voxels:,} ({percent_pasted:.1f}% of input)')
    
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
Adjust the extent of an image by either cropping or expanding.

This tool modifies the spatial extent of an image using either relative offsets or 
absolute extent values. Offsets are applied symmetrically to both sides of each dimension.

OFFSET MODE (--offset X Y Z):
  - Positive values: Expand the image (pad with filler value)
  - Negative values: Crop the image
  - Values are applied symmetrically on both sides of each dimension
  - Negative offsets cannot equal or exceed half the dimension size
  
  Example: --offset 10 10 0 adds 10 voxels to each side in X and Y dimensions
  Example: --offset -5 -5 -5 removes 5 voxels from each side in all dimensions

EXTENT MODE (--extent xmin xmax ymin ymax zmin zmax):
  - Define absolute voxel coordinates for the new extent
  - Can specify both expansion (negative min or max beyond image) and cropping
  - Coordinates are in voxel index space
  
  Example: --extent -10 100 -10 100 0 50 for a 90x90x50 image starting at [-10,-10,0]

The two modes are mutually exclusive - use one or the other, not both.
Areas outside the original image are filled with the filler value (default: 0).
'''
    epilog='''
Example calls: 
ogoImageChangeExtent input.nii.gz output.nii.gz --offset 10 10 10
ogoImageChangeExtent input.nii.gz output.nii.gz --offset -5 -5 0
ogoImageChangeExtent input.nii.gz output.nii.gz --extent -10 200 -10 200 0 100
ogoImageChangeExtent input.nii.gz output.nii.gz -off 20 20 20 --filler -1000
ogoImageChangeExtent input.nii.gz output.nii.gz -ext 0 150 0 150 10 90 -ow -q
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageChangeExtent",
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
    parser.add_argument('--offset', '-off', type=int, nargs=3, default=None, metavar=('X', 'Y', 'Z'),
                        help='Relative offset for each dimension (3 integers); positive=expand, negative=crop')
    parser.add_argument('--extent', '-ext', type=int, nargs=6, default=None, metavar=('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'),
                        help='Absolute extent in voxel coordinates (6 integers)')
    parser.add_argument('--filler', '-f', type=int, default=0, metavar='VALUE',
                        help='Value to use for padding/filling (default: %(default)s)')
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
    if args.input_file is None or args.output_file is None:
        parser.error('input_file and output_file are required')
    
    if not args.quiet:
        ogo.message(echo_arguments('ImageChangeExtent', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = ImageChangeExtent(**run_args)
    return result

if __name__ == '__main__':
    main()
