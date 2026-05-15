# /------------------------------------------------------------------------------+
# | 24-APR-2026                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Post-processing script for ML-based calibration phantom rod segmentation.
# Takes the raw CT scan and the ML model's binary rod segmentation (all rods 
# labeled 1) and automatically identifies the correct calibration phantom by 
# fitting the measured rod HU values against known phantom densities from ogo.
# Rods are then relabeled with their correct phantom-specific labels.

# Imports
import argparse
import inspect
import re

import numpy as np
import nibabel as nib
import scipy as sp
from ogo.calib.standard_calibration import StandardCalibration
from ogo.calib.mindways_calibration import MindwaysCalibration
import ogo.util.Helper as ogo
from ogo.util.echo_arguments import echo_arguments
from ogo.util.cli_tools import (
    check_overwrite,
    validate_input_file
)


#------------------------------------------------------------------------------+
# Helper Functions
def get_known_phantoms():
    """Return a list of all known phantom names."""
    source = inspect.getsource(ogo.get_phantom)
    known_phantoms = re.findall(r"phantom_type in '([^']+)'", source)
    return known_phantoms


def get_all_phantoms():
    """Return a list of all phantom dicts."""
    known_phantoms = get_known_phantoms()
    return [ogo.get_phantom(name) for name in known_phantoms]


def extract_rod_statistics(raw_array, rod_array, quiet):
    """Extract connected components and compute mean HU for each rod."""
    labeled_array, num_features = sp.ndimage.label(rod_array)
    n_rods = num_features
    
    if not quiet:
        ogo.message(f'Found {n_rods} rods')
    
    rod_stats = []
    
    for rod_id in range(1, n_rods + 1):
        rod_mask = (labeled_array == rod_id)
        mean_hu = raw_array[rod_mask].mean()
        
        rod_stats.append({
            "rod_id": rod_id,
            "mean_hu": mean_hu,
            "mask": rod_mask
        })
        
        if not quiet:
            ogo.message(f'  Rod {rod_id}: mean HU = {mean_hu:.2f}')
    
    # Sort rods by mean HU ascending to match to phantom densities later
    rod_stats.sort(key=lambda x: x["mean_hu"])
    
    return rod_stats


def fit_phantom(rod_stats, phantom_dict, quiet):
    """Fit phantom calibration to rod statistics."""
    hu_values = [rod["mean_hu"] for rod in rod_stats]
    densities = phantom_dict["densities"]
    
    if len(hu_values) != len(densities):
        return None
    
    # Sort densities ascending and reorder HU to match
    paired = sorted(zip(densities, hu_values), key=lambda x: x[0])
    densities_sorted = [d for d, h in paired]
    hu_sorted = [h for d, h in paired]
    
    try:
        if phantom_dict["type"] == "k2hpo4":
            water = phantom_dict["h2o_densities"]
            
            # Sort all three lists together by density ascending
            triples = sorted(zip(densities, hu_values, water), key=lambda x: x[0])
            densities_sorted = [d for d, h, w in triples]
            hu_sorted = [h for d, h, w in triples]
            water_sorted = [w for d, h, w in triples]
            
            calibrator = MindwaysCalibration()
            calibrator.fit(hu_sorted, densities_sorted, water_sorted)
        else:
            calibrator = StandardCalibration()
            calibrator.fit(hu_sorted, densities_sorted)
        
        return calibrator
    
    except Exception as e:
        if not quiet:
            ogo.message(f'  Fit failed for {phantom_dict["name"]}: {e}')
        return None


def rank_phantoms(rod_stats, all_phantoms, quiet):
    """Try fitting every known phantom and rank by R^2."""
    n_rods = len(rod_stats)
    results = []
    
    for phantom_dict in all_phantoms:
        if phantom_dict["number_rods"] != n_rods:
            continue
        
        calibrator = fit_phantom(rod_stats, phantom_dict, quiet)
        
        if calibrator is not None:
            r2 = calibrator.r_value ** 2
            if not quiet:
                ogo.message(f'  {phantom_dict["name"]}: R^2 = {r2:.4f}')
            results.append((phantom_dict, calibrator, r2))
    
    results.sort(key=lambda x: x[2], reverse=True)
    return results


def assign_labels(rod_stats, best_phantom, quiet):
    """Assign labels by matching HU rank to density rank."""
    density_label_pairs = sorted(
        zip(best_phantom["densities"], best_phantom["rod_labels"]),
        key=lambda x: x[0]
    )
    
    if not quiet:
        ogo.message('\nLabel assignments:')
    
    for rod, (density, label) in zip(rod_stats, density_label_pairs):
        rod["assigned_label"] = label
        rod["assigned_density"] = density
        if not quiet:
            ogo.message(f'  Rod {rod["rod_id"]}: mean HU = {rod["mean_hu"]:.2f} -> label {label} ({density} mg/cc)')
    
    return rod_stats


def save_labeled_nifti(rod_stats, reference_img, output_file, quiet):
    """Build output label image with correct rod labels and save."""
    output_array = np.zeros(reference_img.shape, dtype=np.int16)
    
    for rod in rod_stats:
        output_array[rod["mask"]] = rod["assigned_label"]
    
    output_img = nib.Nifti1Image(output_array, affine=reference_img.affine, header=reference_img.header)
    nib.save(output_img, output_file)
    
    if not quiet:
        ogo.message(f'\nSaved labeled output to: {output_file}')


#------------------------------------------------------------------------------+
# Main function
def AssignPhantomLabels(raw_nifti_file, rod_nifti_file, output_file, overwrite, quiet):
    """Assign correct labels to segmented phantom rods based on HU values.
    
    This program automatically identifies the calibration phantom from segmented
    phantom rods by fitting measured rod HU values against known phantom densities.
    Rods are then relabeled with their correct phantom-specific labels.
    """
    
    # Check if output exists and should overwrite
    if not check_overwrite(output_file, overwrite):
        return 1
    
    # Load images - raw NIfTI image and segmented image of rods
    if not quiet:
        ogo.message(f'Reading raw image: {raw_nifti_file}')
    
    raw_img = nib.load(raw_nifti_file)
    raw_array = raw_img.get_fdata()
    
    if not quiet:
        ogo.message(f'Reading rod segmentation: {rod_nifti_file}')
    
    rod_img = nib.load(rod_nifti_file)
    rod_array = rod_img.get_fdata()
    
    # Extract rod statistics
    rod_stats = extract_rod_statistics(raw_array, rod_array, quiet)
    
    if len(rod_stats) == 0:
        ogo.message('ERROR: No rods found in segmentation.')
        return 1
    
    # Get all known phantoms and rank by fit
    all_phantoms = get_all_phantoms()
    ranked = rank_phantoms(rod_stats, all_phantoms, quiet)
    
    if not ranked:
        ogo.message(f'ERROR: No known phantom matched {len(rod_stats)} rods. Check segmentation.')
        return 1
    
    # Select best matching phantom
    best_phantom, best_calibrator, best_r2 = ranked[0]
    if not quiet:
        ogo.message(f'\nBest match: {best_phantom["name"]} (R^2 = {best_r2:.4f})')
    
    if best_r2 < 0.99:
        ogo.message(f'  WARNING: R^2 = {best_r2:.4f} is low — check segmentation')
    
    # Assign labels and save output
    rod_stats = assign_labels(rod_stats, best_phantom, quiet)
    save_labeled_nifti(rod_stats, rod_img, output_file, quiet)
    
    if not quiet:
        ogo.message('[DONE]')
    
    return 0


def main():
    """Setup argument parsing and run program."""
    
    description = '''Assign correct labels to segmented calibration phantom rods.

Post-processing script for ML-based calibration phantom rod segmentation.
Takes the raw CT scan and the ML model's binary rod segmentation (all rods 
labeled 1) and automatically identifies the correct calibration phantom by 
fitting the measured rod HU values against known phantom densities. Rods are 
then relabeled with their correct phantom-specific labels.

Output: <rod_segmentation>_labeled.nii.gz — rod mask with correct phantom labels
'''
    
    epilog = '''Example calls:
ogoAssignPhantomLabels raw.nii.gz rods.nii.gz output.nii.gz
ogoAssignPhantomLabels raw.nii.gz rods.nii.gz output.nii.gz --quiet
'''
    
    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoAssignPhantomLabels",
        description=description,
        epilog=epilog
    )
    
    parser.add_argument('raw_nifti_file', type=validate_input_file(['.nii', '.nii.gz']),
                        help='Raw CT image file (*.nii, *.nii.gz)')
    parser.add_argument('rod_nifti_file', type=validate_input_file(['.nii', '.nii.gz']),
                        help='Rod segmentation file (*.nii, *.nii.gz)')
    parser.add_argument('output_file', type=str,
                        help='Output labeled rod segmentation file (*.nii, *.nii.gz)')
    parser.add_argument('--overwrite', '-ow', action='store_true',
                        help='Overwrite output without asking')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress informational output to stdout')
    
    print()
    
    # Parse and display
    args = parser.parse_args()
    
    if not args.quiet:
        ogo.message(echo_arguments('AssignPhantomLabels', vars(args)))
    
    # Run program
    return AssignPhantomLabels(**vars(args))


if __name__ == '__main__':
    exit(main())
