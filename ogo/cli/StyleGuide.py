# /------------------------------------------------------------------------------+
# | 29-MAR-2026                                                                  |
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
    check_label_lengths,
    check_image_dimensions_match,
    validate_integer_range,
    validate_integer_min,
    validate_float_range,
    validate_float_min,
    validate_input_file,
    validate_output_file,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)

#------------------------------------------------------------------------------+
# Main function
def StyleGuide(input_file, mask_file, output_file, from_labels, to_labels, labels_to_keep, collection, threshold, slice_percent, bins, header, overwrite, quiet):

    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Check that from_labels and to_labels have the same length
    if not check_label_lengths(from_labels, to_labels):
        return 1
    
    # Check that at least some labels are defined
    if not from_labels and not to_labels and not labels_to_keep:
        ogo.message('ERROR: No labels defined.')
        ogo.message('       Define --labels_to_keep or --from_labels & --to_labels')
        ogo.message('       There is nothing to do.')
        return 1
    
    # Check if threshold was explicitly set
    if threshold is not None:
        ogo.message(f'  user set threshold:       {threshold}')
    else:
        threshold = 1.0  # Set default value programmatically
    
    # Check if slice_percent was explicitly set
    if slice_percent is not None:
        ogo.message(f'  user set slice_percent:   {slice_percent}')
    else:
        slice_percent = 1.0  # Set default value programmatically
    
    # Check if bins was explicitly set
    if bins is not None:
        ogo.message(f'  user set bins:            {bins}')
    else:
        bins = 10  # Set default value programmatically
    
    # Display label information
    if not quiet:
        ogo.message(f'  from_labels:     {from_labels}')
        ogo.message(f'  to_labels:       {to_labels}')
        ogo.message(f'  labels_to_keep:  {labels_to_keep}')
        ogo.message(f'  threshold:       {threshold}')
        ogo.message(f'  slice_percent:   {slice_percent}')
        ogo.message(f'  bins:            {bins}')
    
    # Read input file based on extension
    if not quiet:
        ogo.message(f'Reading input file: {input_file}')
    
    input_image = sitk.ReadImage(input_file)
    
    if not quiet:
        print_image_info(input_file, input_image)
    
    # Read mask file based on extension
    if not quiet:
        ogo.message(f'Reading mask file: {mask_file}')
    
    mask_image = sitk.ReadImage(mask_file, sitk.sitkUInt8)

    if not quiet:
        print_image_info(mask_file, mask_image)
    
    # Check that images have matching dimensions
    if not check_image_dimensions_match(input_image, mask_image, "input_file", "mask_file"):
        return 1
    
    # Function logic goes here
        
    if not quiet:
        ogo.message('[DONE]')
    
    return 0

def main():
    # Setup description
    description='''
This program demonstrates the coding styles and best practices for Ogo CLI tools.

Basically it does nothing... but it shows how to set up a CLI with argparse, 
how to validate inputs, and how to structure the code for readability and 
maintainability.

STYLE GUIDE RULES:
  1. Use the TYPE option in argparse to check valid inputs whenever possible.
  2. Minimize use of positional arguments for mandatory inputs only.
  3. Use flags with long names, but include short names in addition if needed (see below).
  4. Return 0 for successful completion, return 1 for early exit.
  5. Formatting:
       – Use 'snake_case' for variables and function names (e.g., input_file).
       – Use 'PascalCase' for class names (e.g., ImageProcessor).
       - Use 'UPPER_SNAKE_CASE' for constants (e.g., MAX_ITERATIONS).
       - Do NOT use 'camelCase' (e.g., inputLabels, outputLabels).

EXAMPLE VARIABLE NAMES:
  input_file / output_file / mask_file / calib_file
  input_path / output_path
  from_labels / to_labels / labels_to_keep
  overwrite

COMMON ARGUMENT OPTIONS:
  -h,  --help      Shows user information about CLI
  -ow, --overwrite Overwrite output if it exists
  -o,  --output    Output file (e.g., sort, gcc)
  -q,  --quiet     Display less output (useful for scripts)
  -d,  --debug     Show debugging output (verbose)
  -f,  --force     Force destructive actions without confirmation
       --json      Display JSON output
  -n,  --dry_run   Describe changes without executing
  -a,  --all       All (e.g., ps, fetchmail)
  -u,  --user      User (e.g., ps, ssh)
  -v,  --version   Version

We can use JSON files instead of command line arguments for more complex configurations.
  - To output current arguments as JSON: --json
  - To load arguments from JSON file: --json my_config.json  
'''
    epilog='''
Example calls: 
ogoStyleGuide input.nii.gz output.nii.gz -fl 1 2 3 -tl 10 20 30
ogoStyleGuide input.nii.gz output.nii.gz -lk 5 10 15 --quiet
ogoStyleGuide input.nii.gz output.nii.gz -lk 12 3 --json > my_config.json
ogoStyleGuide --json my_config.json
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoStyleGuide",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('mask_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']), 
                        help='Mask image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.csv', '.txt']), 
                        help='Output file (*.csv, *.txt)')
    parser.add_argument('--from_labels','-fl', type=validate_integer_range(0, 255), nargs='*', default=[], metavar='ID', help='Target label input IDs; space separated (e.g. 1 2 3)')
    parser.add_argument('--to_labels','-tl', type=validate_integer_range(0, 255), nargs='*', default=[], metavar='ID', help='Target label output IDs; space separated (e.g. 4 5 6)')
    parser.add_argument('--labels_to_keep','-lk', type=validate_integer_range(0, 255), nargs='*', default=[], metavar='ID', help='Labels to keep (others removed); (e.g. 7 8 9)')
    parser.add_argument('--collection', default='all', choices=['all', 'ossai', 'skeleton', 'cardio', 'organs', 'muscle', 'spine', 'gastro'],
                                                           help='Select collection (default: %(default)s)')
    parser.add_argument('--threshold', '-t', type=validate_float_min(1.0), default=None, help='Threshold value (must be >= 10.0, default: 1.0)')
    parser.add_argument('--slice_percent', '-s', type=validate_float_range(1.0, 5.0), default=None, help='Slice percentage (must be between 1.0 and 5.0, default: 1.0)')
    parser.add_argument('--bins', '-b', type=validate_integer_min(1), default=None, help='Number of bins (must be >= 1, default: 10)')
    parser.add_argument('--header','-H', action='store_true', help='Include header in output')
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
        ogo.message(echo_arguments('StyleGuide', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    StyleGuide(**run_args)

if __name__ == '__main__':
    main()
