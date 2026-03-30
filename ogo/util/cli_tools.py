# /------------------------------------------------------------------------------+
# | 29-MAR-2026                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

"""
CLI tools and validation utility functions for argument parsing and file operations.
"""

# Imports
import os
import math
import json
import argparse
import SimpleITK as sitk

def check_overwrite(output_file, overwrite):
    """
    Check if output file can be overwritten.
    
    Args:
        output_file: Path to the output file
        overwrite: Whether to overwrite without asking
    
    Returns:
        True if it's safe to proceed (file doesn't exist, overwrite is True, or user confirmed)
        False if user declined to overwrite
    """
    if not overwrite and os.path.isfile(output_file):
        result = input(f'File "{output_file}" already exists. Overwrite? [y/n]: ')
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            return False
    return True

def check_label_lengths(from_labels, to_labels):
    """
    Check that from_labels and to_labels have the same length.
    
    Args:
        from_labels: List of source labels
        to_labels: List of target labels
    
    Returns:
        True if lengths match or both are empty
        False if lengths don't match (prints error and exits)
    """
    if len(from_labels) != len(to_labels):
        print(f'ERROR: Number of from_labels must equal number of to_labels.')
        print(f'       [from_labels #{len(from_labels)} != to_labels #{len(to_labels)}]')
        return False
    return True

def check_image_dimensions_match(image1, image2, name1="image1", name2="image2"):
    """
    Check if two images have the same dimensions.
    
    Args:
        image1: First SimpleITK image
        image2: Second SimpleITK image
        name1: Name/description of first image for error message
        name2: Name/description of second image for error message
    
    Returns:
        True if dimensions match, False otherwise
    """
    size1 = image1.GetSize()
    size2 = image2.GetSize()
    
    if size1 != size2:
        print(f'ERROR: Image dimensions do not match.')
        print(f'       {name1}: {size1}')
        print(f'       {name2}: {size2}')
        return False
    return True

def validate_integer_range(min_val, max_val):
    """
    Factory function that returns a validator for integers within a range.
    
    Args:
        min_val: Minimum allowed value (inclusive)
        max_val: Maximum allowed value (inclusive)
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(value):
        """Validates that the integer is within the specified range."""
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'Invalid integer value: "{value}"')
        
        if ivalue < min_val or ivalue > max_val:
            raise argparse.ArgumentTypeError(
                f'Value must be between {min_val} and {max_val}, got {ivalue}'
            )
        return ivalue
    
    return validator

def validate_integer_min(min_val):
    """
    Factory function that returns a validator for integers greater than or equal to a minimum value.
    
    Args:
        min_val: Minimum allowed value (inclusive)
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(value):
        """Validates that the integer is greater than or equal to the minimum value."""
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'Invalid integer value: "{value}"')
        
        if ivalue < min_val:
            raise argparse.ArgumentTypeError(
                f'Value must be greater than or equal to {min_val}, got {ivalue}'
            )
        return ivalue
    
    return validator

def validate_float_range(min_val, max_val):
    """
    Factory function that returns a validator for floats within a range.
    
    Args:
        min_val: Minimum allowed value (inclusive)
        max_val: Maximum allowed value (inclusive)
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(value):
        """Validates that the float is within the specified range."""
        try:
            fvalue = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'Invalid float value: "{value}"')
        
        if fvalue < min_val or fvalue > max_val:
            raise argparse.ArgumentTypeError(
                f'Value must be between {min_val} and {max_val}, got {fvalue}'
            )
        return fvalue
    
    return validator

def validate_path_exists(value):
    """
    Validator function to check if a path or its parent directory exists.
    
    This function validates two cases:
    1. The path itself exists (file or directory)
    2. The parent directory exists (useful for output files that don't exist yet)
    
    Args:
        value: Path string to validate
    
    Returns:
        The path if it or its parent directory exists
        
    Raises:
        argparse.ArgumentTypeError: If neither the path nor its parent directory exists
    """
    # Check if the path exists as-is
    if os.path.exists(value):
        return value
    
    # If not, check if the parent directory exists (for output files)
    parent_dir = os.path.dirname(value)
    
    # Handle case where dirname returns empty string (current directory)
    if parent_dir == '':
        parent_dir = '.'
    
    if os.path.exists(parent_dir):
        return value
    
    # Neither path nor parent directory exists
    raise argparse.ArgumentTypeError(
        f'Path does not exist: "{value}". Parent directory "{parent_dir}" also does not exist.'
    )

def validate_input_file(allowed_extensions):
    """
    Factory function that returns a validator for input files.
    Checks both file existence and valid extensions.
    
    Args:
        allowed_extensions: List or tuple of allowed file extensions (e.g., ['.nii', '.nii.gz', '.aim'])
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(filepath):
        """Validates input file existence and extension."""
        # Check if file exists
        if not os.path.isfile(filepath):
            raise argparse.ArgumentTypeError(
                f'File does not exist: "{filepath}"'
            )
        
        # Check file extension
        filepath_lower = filepath.lower()
        if not any(filepath_lower.endswith(ext) for ext in allowed_extensions):
            # Get the actual extension
            if '.' in filepath:
                actual_ext = '.' + filepath.rsplit('.', 1)[1]
                if filepath_lower.endswith('.nii.gz'):
                    actual_ext = '.nii.gz'
            else:
                actual_ext = 'no extension'
            
            raise argparse.ArgumentTypeError(
                f'Invalid file extension: "{actual_ext}". '
                f'Allowed extensions: {", ".join(allowed_extensions)}'
            )
        
        return filepath
    
    return validator

def validate_output_file(allowed_extensions):
    """
    Factory function that returns a validator for output files.
    Checks valid extensions only (does not check file existence).
    
    Args:
        allowed_extensions: List or tuple of allowed file extensions (e.g., ['.nii', '.nii.gz', '.csv', '.txt'])
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(filepath):
        """Validates output file extension."""
        # Check file extension
        filepath_lower = filepath.lower()
        if not any(filepath_lower.endswith(ext) for ext in allowed_extensions):
            # Get the actual extension
            if '.' in filepath:
                actual_ext = '.' + filepath.rsplit('.', 1)[1]
                if filepath_lower.endswith('.nii.gz'):
                    actual_ext = '.nii.gz'
            else:
                actual_ext = 'no extension'
            
            raise argparse.ArgumentTypeError(
                f'Invalid output file extension: "{actual_ext}". '
                f'Allowed extensions: {", ".join(allowed_extensions)}'
            )
        
        return filepath
    
    return validator

def validate_float_min(min_val):
    """
    Factory function that returns a validator for floats greater than or equal to a minimum value.
    
    Args:
        min_val: Minimum allowed value (inclusive)
    
    Returns:
        A validator function for use with argparse type parameter
    """
    def validator(value):
        """Validates that the float is greater than or equal to the minimum value."""
        try:
            fvalue = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f'Invalid float value: "{value}"')
        
        if fvalue < min_val:
            raise argparse.ArgumentTypeError(
                f'Value must be greater than or equal to {min_val}, got {fvalue}'
            )
        return fvalue
    
    return validator

def load_config_from_json(json_file):
    """
    Load configuration from JSON file.
    
    Args:
        json_file: Path to JSON configuration file
    
    Returns:
        Dictionary of configuration parameters
    """
    with open(json_file, 'r') as f:
        return json.load(f)

def save_config_to_json(args, json_file):
    """
    Save configuration to JSON file or output to stdout.
    
    Args:
        args: Parsed arguments namespace
        json_file: Path to output JSON file, or None to output to stdout
    """
    config = vars(args).copy()
    # Remove json-related arguments
    config.pop('json', None)
    
    if json_file is None:
        # Output to stdout
        print(json.dumps(config, indent=2))
    else:
        # Validate file extension
        if not json_file.endswith('.json'):
            raise ValueError(f'Configuration file must have .json extension: {json_file}')
        
        # Save to file
        with open(json_file, 'w') as f:
            json.dump(config, f, indent=2)
        print(f'Configuration saved to: {json_file}')

def print_image_info(infile, image):
    """
    Print formatted image information including dimensions, spacing, and statistics.
    
    Args:
        infile: Path to the image file
        image: SimpleITK image object
    """
    guard = '!-------------------------------------------------------------------------------'
    phys_dim = [x * y for x, y in zip(image.GetSize(), image.GetSpacing())]
    position = [math.floor(x / y) for x, y in zip(image.GetOrigin(), image.GetSpacing())]
    size = os.path.getsize(infile)
    names = ['Bytes', 'KBytes', 'MBytes', 'GBytes']
    n_image_voxels = image.GetSize()[0] * image.GetSize()[1] * image.GetSize()[2]
    voxel_volume = image.GetSpacing()[0] * image.GetSpacing()[1] * image.GetSpacing()[2]
    
    # Calculate file size with proper unit
    i = 0
    while size > 1024 and i < len(names) - 1:
        i += 1
        size = size / 1024.0

    # Get image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    
    # Print formatted information
    print(guard)
    print('!>')
    print('!> File                       {}'.format(os.path.basename(infile)))
    print('!> dim                            {:>8}  {:>8}  {:>8}'.format(*image.GetSize()))
    print('!> off                            {:>8}  {:>8}  {:>8}'.format('-', '-', '-'))
    print('!> pos                            {:>8}  {:>8}  {:>8}'.format(*position))
    print('!> element size in mm             {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*image.GetSpacing()))
    print('!> phys dim in mm                 {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*phys_dim))
    print('!>')
    print('!> Number of voxels           {:>12,}'.format(n_image_voxels))
    print('!> Voxel volume (mm^3)        {:>12.6f}'.format(voxel_volume))
    print('!> Type of data               {}'.format(image.GetPixelIDTypeAsString()))
    print('!> Image min/max              {:>12.2f} / {:.2f}'.format(stats.GetMinimum(), stats.GetMaximum()))
    print('!> Total memory size          {:.1f} {}'.format(size, names[i]))
    print(guard)
