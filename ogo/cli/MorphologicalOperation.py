# /------------------------------------------------------------------------------+
# | 09-MAR-2023                                                                  |
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
    validate_integer_range,
    validate_integer_min,
    load_config_from_json,
    save_config_to_json,
    print_image_info
)


def MorphologicalOperation(input_file, output_file, operation_labels, dilate, erode, morphological, overwrite, quiet):
    
    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Validate dilate and erode parameters
    if not dilate:
        dilate = []
    if not erode:
        erode = []
    
    # Check that at least one operation is specified
    if not dilate and not erode:
        ogo.message('[ERROR] At least one of --dilate or --erode must be specified.')
        return 1
    
    # Validate all dilate values are positive
    for d in dilate:
        if d <= 0:
            ogo.message(f'[ERROR] All dilate kernel values must be positive. Found: {d}')
            return 1
    
    # Validate all erode values are positive
    for e in erode:
        if e <= 0:
            ogo.message(f'[ERROR] All erode kernel values must be positive. Found: {e}')
            return 1
    
    # Read input
    if not quiet:
        ogo.message(f'Reading input image:')
        ogo.message(f'  {input_file}')
    
    input_image = sitk.ReadImage(input_file, sitk.sitkUInt8)
    
    if not quiet:
        print_image_info(input_file, input_image, brief=True)
    
    # Get all available labels in the input image
    label_filter = sitk.LabelShapeStatisticsImageFilter()
    label_filter.Execute(input_image)
    available_labels = list(label_filter.GetLabels())
    
    if not quiet:
        ogo.message(f'Available labels in input image: {available_labels}')
    
    # Set the labels on which we will perform morphological operations
    if not operation_labels:
        operation_labels = available_labels
        if not quiet:
            ogo.message('No operation_labels specified - will process all labels')
    else:
        # Filter operation_labels to only include labels that exist in the image
        original_operation_labels = operation_labels.copy()
        valid_operation_labels = []
        invalid_labels = []
        
        # Track which indices have valid labels (for per-label kernels)
        valid_indices = []
        
        for idx, label in enumerate(original_operation_labels):
            if label in available_labels:
                valid_operation_labels.append(label)
                valid_indices.append(idx)
            else:
                invalid_labels.append(label)
        
        # Warn about invalid labels
        if invalid_labels:
            if not quiet:
                ogo.message(f'[WARNING] Label(s) {invalid_labels} not found in input image and will be skipped.')
                ogo.message(f'          Available labels: {available_labels}')
        
        # Update operation_labels to only valid labels
        operation_labels = valid_operation_labels
        
        # If no valid labels remain, nothing to do
        if not operation_labels:
            ogo.message('[ERROR] None of the specified operation_labels exist in the input image.')
            ogo.message(f'        Available labels: {available_labels}')
            return 1
        
        # Filter dilate and erode to match valid labels if they were specified per-label
        if len(dilate) == len(original_operation_labels):
            dilate = [dilate[i] for i in valid_indices]
        
        if len(erode) == len(original_operation_labels):
            erode = [erode[i] for i in valid_indices]
    
    # Expand dilate and erode to match number of labels
    num_labels = len(operation_labels)
    
    if len(dilate) == 1:
        dilate_kernels = dilate * num_labels
    elif len(dilate) == num_labels:
        dilate_kernels = dilate
    elif len(dilate) == 0:
        dilate_kernels = [0] * num_labels
    else:
        ogo.message('[ERROR] Number of --dilate values must be either 1 or match number of operation labels.')
        ogo.message(f'        operation_labels: {num_labels}, dilate values: {len(dilate)}')
        return 1
    
    if len(erode) == 1:
        erode_kernels = erode * num_labels
    elif len(erode) == num_labels:
        erode_kernels = erode
    elif len(erode) == 0:
        erode_kernels = [0] * num_labels
    else:
        ogo.message('[ERROR] Number of --erode values must be either 1 or match number of operation labels.')
        ogo.message(f'        operation_labels: {num_labels}, erode values: {len(erode)}')
        return 1
    
    # Warn if both dilate and erode are specified without explicit morphological parameter
    if dilate and erode and morphological is None:
        if not quiet:
            ogo.message('[WARNING] Both --dilate and --erode specified. Dilation will be applied first, then erosion.')
            ogo.message('          This is equivalent to a morphological closing operation.')
            ogo.message('          To explicitly control operation order, use --morphological opening or --morphological closing.')
    
    # Warn if morphological operation selected but parameters not aligned
    if morphological == 'opening' and not erode:
        if not quiet:
            ogo.message('[WARNING] Morphological opening selected but --erode not specified.')
            ogo.message('          Only dilation will be applied.')
    
    if morphological == 'closing' and not dilate:
        if not quiet:
            ogo.message('[WARNING] Morphological closing selected but --dilate not specified.')
            ogo.message('          Only erosion will be applied.')
    
    # Display operation info
    if not quiet:
        if morphological == 'opening':
            ogo.message('Morphological operation: opening (erode then dilate)')
        elif morphological == 'closing':
            ogo.message('Morphological operation: closing (dilate then erode)')
        else:
            ogo.message('Morphological operation: dilate then erode (default)')
    
    # Initialize output as empty segmentation
    output_seg = input_image < 0  # All zeros
    
    # Process each label
    for label in available_labels:
        # Create binary mask for this label
        label_mask = input_image == label
        
        if label in operation_labels:
            # Find the kernels for this label
            label_idx = operation_labels.index(label)
            current_dilate = dilate_kernels[label_idx]
            current_erode = erode_kernels[label_idx]
            
            # Apply morphological operations based on morphological parameter
            processed_mask = label_mask
            
            if morphological == 'opening':
                # Opening: erode then dilate
                if current_erode > 0:
                    erode_filter = sitk.BinaryErodeImageFilter()
                    erode_filter.SetKernelRadius(current_erode)
                    processed_mask = erode_filter.Execute(processed_mask)
                
                if current_dilate > 0:
                    dilate_filter = sitk.BinaryDilateImageFilter()
                    dilate_filter.SetKernelRadius(current_dilate)
                    processed_mask = dilate_filter.Execute(processed_mask)
                
                if not quiet:
                    ogo.message(f'  Label {label:3d}: opening (erode={current_erode}, dilate={current_dilate})')
            
            elif morphological == 'closing':
                # Closing: dilate then erode
                if current_dilate > 0:
                    dilate_filter = sitk.BinaryDilateImageFilter()
                    dilate_filter.SetKernelRadius(current_dilate)
                    processed_mask = dilate_filter.Execute(processed_mask)
                
                if current_erode > 0:
                    erode_filter = sitk.BinaryErodeImageFilter()
                    erode_filter.SetKernelRadius(current_erode)
                    processed_mask = erode_filter.Execute(processed_mask)
                
                if not quiet:
                    ogo.message(f'  Label {label:3d}: closing (dilate={current_dilate}, erode={current_erode})')
            
            else:
                # Default: dilate then erode
                if current_dilate > 0:
                    dilate_filter = sitk.BinaryDilateImageFilter()
                    dilate_filter.SetKernelRadius(current_dilate)
                    processed_mask = dilate_filter.Execute(processed_mask)
                
                if current_erode > 0:
                    erode_filter = sitk.BinaryErodeImageFilter()
                    erode_filter.SetKernelRadius(current_erode)
                    processed_mask = erode_filter.Execute(processed_mask)
                
                if not quiet:
                    ogo.message(f'  Label {label:3d}: dilate={current_dilate}, erode={current_erode}')
            
            morphed_mask = processed_mask
        else:
            # No operation for this label - keep as is
            morphed_mask = label_mask
            if not quiet:
                ogo.message(f'  Label {label:3d}: unchanged')
        
        # Handle overlaps - existing labels have priority
        existing_seg = output_seg > 0
        new_label_mask = morphed_mask > 0
        
        # Cast to integers and find overlap
        existing_seg_int = sitk.Cast(existing_seg, sitk.sitkUInt8)
        new_label_mask_int = sitk.Cast(new_label_mask, sitk.sitkUInt8)
        overlap = (existing_seg_int + new_label_mask_int) == 2
        
        # Mask out overlapping regions from new label
        no_overlap_mask = sitk.Cast(overlap == 0, sitk.sitkUInt8)
        this_label_clean = sitk.Mask(morphed_mask, no_overlap_mask)
        
        # Add this label to the output
        output_seg = output_seg + sitk.Cast(label * (this_label_clean > 0), sitk.sitkUInt8)
    
    # Write output
    if not quiet:
        ogo.message(f'Writing output to: {output_file}')
    
    sitk.WriteImage(output_seg, output_file)
    
    if not quiet:
        print_image_info(output_file, output_seg, brief=True)
        ogo.message('[DONE]')
    
    return 0


def main():
    description = '''
Perform morphological operations on labeled images.

This tool applies binary morphological operations (dilation and erosion) to specific 
labels in a segmented image. Each label is processed independently, and overlaps are 
resolved by giving priority to labels processed first (in the order they appear in the image).

MORPHOLOGICAL OPERATIONS:
  - opening:  Erode then dilate (removes small protrusions, smooths from outside)
  - closing:  Dilate then erode (fills small holes, smooths from inside)
  - Default:  Dilate then erode (if both --dilate and --erode specified)

DILATE AND ERODE:
  Specify kernel radii for dilation and erosion operations:
  - Single value: applies to all processed labels
  - Multiple values: must match number of operation_labels (one kernel per label)
  - Kernel radius determines structuring element size (1=3x3x3, 2=5x5x5, etc.)

LABEL PROCESSING:
  - If no operation_labels specified, all labels in the image are processed
  - Labels are processed in the order they appear in the image
  - Overlaps are resolved by keeping the first label (earlier labels have priority)

JSON CONFIGURATION:
  - To output current arguments as JSON: --json
  - To load arguments from JSON file: --json my_config.json
'''

    epilog = '''
Example calls: 
ogoMorphologicalOperation input.nii.gz output.nii.gz --dilate 2
ogoMorphologicalOperation input.nii.gz output.nii.gz --erode 3 --operation_labels 1 2 3
ogoMorphologicalOperation input.nii.gz output.nii.gz --dilate 2 --erode 1 --morphological closing
ogoMorphologicalOperation input.nii.gz output.nii.gz --dilate 2 3 4 --erode 1 1 1 --operation_labels 1 2 3
ogoMorphologicalOperation input.nii.gz output.nii.gz --dilate 2 --morphological opening --json > config.json
ogoMorphologicalOperation --json config.json
'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoMorphologicalOperation",
        description=description,
        epilog=epilog
    )
    
    # JSON configuration argument
    parser.add_argument('--json', nargs='?', const='', metavar='FILE',
                        help='Output arguments as JSON to stdout (no file), or load arguments from JSON file')
    
    parser.add_argument('input_file', nargs='?', type=validate_input_file(['.nii', '.nii.gz']),
                        help='Input label image file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', nargs='?', type=validate_output_file(['.nii', '.nii.gz']),
                        help='Output label image file (*.nii, *.nii.gz)')
    parser.add_argument('--operation_labels','-l', type=validate_integer_range(0, 255), nargs='*', default=[], 
                        metavar='LABEL', help='List of labels to process (default: all labels in image)')
    parser.add_argument('--dilate', '-d', type=validate_integer_min(1), nargs='*', default=[], metavar='KERNEL',
                        help='Kernel radius/radii for dilation (single value or one per label)')
    parser.add_argument('--erode', '-e', type=validate_integer_min(1), nargs='*', default=[], metavar='KERNEL',
                        help='Kernel radius/radii for erosion (single value or one per label)')
    parser.add_argument('--morphological', '-mo', choices=['opening', 'closing'], default=None, metavar='TYPE',
                        help='Morphological operation type: opening (erode then dilate) or closing (dilate then erode)')
    parser.add_argument('--overwrite', '-ow', action='store_true', 
                        help='Overwrite output without asking')
    parser.add_argument('--quiet', '-q', action='store_true', 
                        help='Suppress informational output to stdout')

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
        ogo.message(echo_arguments('MorphologicalOperation', vars(args)))

    # Run program (exclude json argument)
    run_args = vars(args).copy()
    run_args.pop('json', None)
    result = MorphologicalOperation(**run_args)
    return result


if __name__ == '__main__':
    main()