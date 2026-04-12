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

def validate_extent(extent_list):
    """
    Validator for image extent specified as [xmin, xmax, ymin, ymax, zmin, zmax].
    Ensures that min < max for each dimension.
    
    Args:
        extent_list: List of 6 integers representing the extent
    
    Returns:
        The validated extent list
        
    Raises:
        argparse.ArgumentTypeError: If validation fails
    """
    if len(extent_list) != 6:
        raise argparse.ArgumentTypeError('set_extent requires exactly 6 integers')
    
    xmin, xmax, ymin, ymax, zmin, zmax = extent_list
    
    if xmin >= xmax:
        raise argparse.ArgumentTypeError(f'xmin ({xmin}) must be less than xmax ({xmax})')
    if ymin >= ymax:
        raise argparse.ArgumentTypeError(f'ymin ({ymin}) must be less than ymax ({ymax})')
    if zmin >= zmax:
        raise argparse.ArgumentTypeError(f'zmin ({zmin}) must be less than zmax ({zmax})')
    
    return extent_list

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
    Print formatted image information including dimensions, spacing, orientation, and statistics.
    
    Args:
        infile: Path to the image file
        image: SimpleITK image object
    """
    import numpy as np
    
    guard = '!-------------------------------------------------------------------------------'
    
    # Get basic image properties
    size = image.GetSize()
    spacing = image.GetSpacing()
    origin = image.GetOrigin()
    direction = image.GetDirection()
    
    # Calculate derived properties
    phys_dim = [x * y for x, y in zip(size, spacing)]
    position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    n_image_voxels = size[0] * size[1] * size[2]
    voxel_volume = spacing[0] * spacing[1] * spacing[2]
    total_physical_volume = phys_dim[0] * phys_dim[1] * phys_dim[2]
    
    # Calculate physical extent (bounding box)
    phys_min = origin
    phys_max = image.TransformContinuousIndexToPhysicalPoint([float(size[0]), float(size[1]), float(size[2])])
    
    # Calculate direction matrix determinant (handedness)
    direction_matrix = np.array(direction).reshape(3, 3)
    determinant = np.linalg.det(direction_matrix)
    handedness = "right-handed" if determinant > 0 else "left-handed"
    
    # Get file size with proper unit
    try:
        file_size = os.path.getsize(infile)
        size_names = ['Bytes', 'KBytes', 'MBytes', 'GBytes']
        size_unit_index = 0
        while file_size > 1024 and size_unit_index < len(size_names) - 1:
            size_unit_index += 1
            file_size = file_size / 1024.0
    except OSError:
        file_size = 0
        size_unit_index = 0
        size_names = ['Bytes', 'KBytes', 'MBytes', 'GBytes']

    # Get image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    
    # Print formatted information
    print(guard)
    print('!>')
    print(f'!> File                       {os.path.basename(infile)}')
    print(f'!> dim                            {size[0]:>8}  {size[1]:>8}  {size[2]:>8}')
    print(f'!> off                            {"-":>8}  {"-":>8}  {"-":>8}')
    print(f'!> pos                            {position[0]:>8}  {position[1]:>8}  {position[2]:>8}')
    print(f'!> element size in mm             {spacing[0]:>8.4f}  {spacing[1]:>8.4f}  {spacing[2]:>8.4f}')
    print(f'!> phys dim in mm                 {phys_dim[0]:>8.4f}  {phys_dim[1]:>8.4f}  {phys_dim[2]:>8.4f}')
    print(f'!> origin in mm                   {origin[0]:>8.4f}  {origin[1]:>8.4f}  {origin[2]:>8.4f}')
    print('!>')
    print(f'!> Physical extent (mm):')
    print(f'!>   min                          {phys_min[0]:>8.4f}  {phys_min[1]:>8.4f}  {phys_min[2]:>8.4f}')
    print(f'!>   max                          {phys_max[0]:>8.4f}  {phys_max[1]:>8.4f}  {phys_max[2]:>8.4f}')
    print('!>')
    print(f'!> Direction matrix ({handedness}):')
    print(f'!>   [{direction[0]:>7.4f}  {direction[1]:>7.4f}  {direction[2]:>7.4f}]')
    print(f'!>   [{direction[3]:>7.4f}  {direction[4]:>7.4f}  {direction[5]:>7.4f}]')
    print(f'!>   [{direction[6]:>7.4f}  {direction[7]:>7.4f}  {direction[8]:>7.4f}]')
    print('!>')
    print(f'!> Number of voxels           {n_image_voxels:>12,}')
    print(f'!> Voxel volume (mm^3)        {voxel_volume:>12.6f}')
    print(f'!> Total physical volume (mm^3) {total_physical_volume:>10.2f}')
    print(f'!> Total physical volume (cm^3) {total_physical_volume/1000.0:>10.4f}')
    print(f'!> Type of data               {image.GetPixelIDTypeAsString()}')
    print(f'!> Image min/max              {stats.GetMinimum():>12.2f} / {stats.GetMaximum():.2f}')
    print(f'!> Total memory size          {file_size:.1f} {size_names[size_unit_index]}')
    print(guard)
