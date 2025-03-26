import os
import cv2
import numpy as np
import SimpleITK as sitk

from collections import defaultdict
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from scipy.ndimage import (
    binary_erosion, binary_dilation, gaussian_filter1d, binary_fill_holes
)

from skimage.filters import gaussian
from skimage.morphology import binary_opening, remove_small_objects
from skimage.measure import label, regionprops
from skimage.feature import peak_local_max, canny
from skimage.draw import disk
from skimage.transform import hough_circle, hough_circle_peaks
import argparse 


def model_labels():
    
    model_label_map = {
        "Mindways Model 3 CT": [111, 112, 113, 114, 115],
        "B-MAS200": [116, 117, 118, 119, 120],
        "Mindways Model 3 QA": [121, 122, 123, 124],
        "QRM-BDC 3-rod": [125, 126, 127],
        "QRM-BDC 6-rod": [128, 129, 130, 131, 132, 133],
        "Image Analysis QCT-3D Plus": [134, 135, 136],
    }
    return model_label_map

def shift_image_by_background_peak(im_arr, mask_arr, hist_bins=2048, z_cutoff=3.0, debug=False):
    """
    Detects the dominant intensity peak in the masked region and shifts the image so this peak is centered at 0.

    Parameters:
        im_arr (ndarray): 3D image array.
        mask_arr (ndarray): 3D binary mask (1 = region of interest, 0 = ignore).
        hist_bins (int): Number of bins for histogram.
        z_cutoff (float): Standard deviation multiplier to define near-zero region (used in optional plotting).
        debug (bool): If True, show histogram and peak fit.

    Returns:
        shifted (ndarray): Image with main peak shifted to 0.
        abs_shifted (ndarray): Absolute value of shifted image.
        mean_est (float): Estimated mean of background peak.
        std_est (float): Estimated std deviation of background peak.
    """
    # Smooth and erode mask to avoid edges
    smoothed = gaussian(im_arr, sigma=1.0, preserve_range=True)
    eroded_mask = binary_erosion(mask_arr)
    masked_arr = smoothed.copy()
    masked_arr[eroded_mask == 0] = np.nan

    # Histogram of valid values
    valid_values = masked_arr[~np.isnan(masked_arr)].ravel()
    hist, bin_edges = np.histogram(valid_values, bins=hist_bins, range=(np.nanmin(valid_values), np.nanmax(valid_values)))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smoothed_hist = gaussian_filter1d(hist, sigma=3)

    # Detect peaks and estimate main peak location
    peak_indices = peak_local_max(smoothed_hist, num_peaks=7, exclude_border=False, threshold_rel=0.01)
    peak_values = bin_centers[peak_indices[:, 0]]
    sorted_peaks = sorted(zip(peak_values, smoothed_hist[peak_indices[:, 0]]), key=lambda x: x[1], reverse=True)
    background_peak = sorted_peaks[0][0]
    background_idx = np.argmin(np.abs(bin_centers - background_peak))

    # Fit Gaussian-like window around peak
    window_size = 50
    start_idx = max(0, background_idx - window_size)
    end_idx = min(len(bin_centers), background_idx + window_size)
    x_window = bin_centers[start_idx:end_idx]
    y_window = smoothed_hist[start_idx:end_idx]
    mean_est = np.average(x_window, weights=y_window)
    std_est = np.sqrt(np.average((x_window - mean_est) ** 2, weights=y_window))

    # Shift image
    shifted = smoothed - mean_est
    abs_shifted = np.abs(shifted)

    return abs_shifted


def hough_detection(image_arr, mask_arr, radius_range=(6, 10), debug=True, nrods=5, param2_start=30, param2_min=5):
    # Shift image so background peak is centered at 0
    shifted  = shift_image_by_background_peak(image_arr, mask_arr)
    
    blurrarr = gaussian(binary_erosion(mask_arr).astype(float), sigma=2.0, preserve_range=True)
    shifted*=blurrarr
    shifted[shifted<50]=0
    
    output_mask = np.zeros_like(image_arr, dtype=np.uint8)
    z_indices = np.where(mask_arr.any(axis=(1, 2)))[0]
    min_radius, max_radius = radius_range
    params = []
    for z in z_indices[5:-5]:
        slice_img = shifted[z].copy()
        slice_mask = mask_arr[z]
        slice_img[slice_mask == 0] = 0

        # Normalize to 8-bit grayscale for OpenCV
        norm_slice = cv2.normalize(slice_img, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        blurred = cv2.GaussianBlur(norm_slice, (5, 5), 1)

        # Try varying param2 from high to low sensitivity
        found_circles = None
        for param2 in range(param2_start, param2_min - 1, -1):
            circles = cv2.HoughCircles(
                blurred,
                cv2.HOUGH_GRADIENT,
                dp=1.2,
                minDist=20,
                param1=30,
                param2=param2,
                minRadius=min_radius,
                maxRadius=max_radius
            )

            if circles is not None and len(circles[0]) >= nrods:
                found_circles = np.uint16(np.around(circles[0][:nrods]))
                break  # exit loop when enough circles are found
        
        params.append(param2)

        # Draw the selected circles
        if found_circles is not None:
            for center_x, center_y, radius in found_circles:
                rr, cc = disk((center_y, center_x), radius, shape=slice_img.shape)
                output_mask[z, rr, cc] = 1

    
    print(f'Using Param2: {np.mean(param2)}')

    return output_mask


def extract_slicewise_centroids(mask):
    """
    Returns a dictionary {z: [(y, x)]} of centroids per slice.
    """
    Z = mask.shape[0]
    centroids_per_slice = defaultdict(list)

    for z in range(Z):
        labeled_slice = label(mask[z])
        props = regionprops(labeled_slice)
        for p in props:
            y, x = p.centroid
            centroids_per_slice[z].append((x, y))  # note: x, y ordering!

    return centroids_per_slice

def estimate_num_components(model):

    model_label_map = model_labels()
    return len(model_label_map[model])

def filter_full_slices(centroids_per_slice, n_expected):
    return {z: centroids for z, centroids in centroids_per_slice.items() if len(centroids) == n_expected}


def cluster_centroids(centroids_per_slice_full, n_components):
    """
    Returns dict: {component_id: [(z, x, y), ...]} with consistent IDs.
    """
    all_points = []
    for z, centroids in centroids_per_slice_full.items():
        for x, y in centroids:
            all_points.append([z, x, y])

    all_points = np.array(all_points)

    # cluster based on x/y (since z varies by slice)
    kmeans = KMeans(n_clusters=n_components, random_state=42).fit(all_points[:, 1:3])
    labels = kmeans.labels_

    component_dict = defaultdict(list)
    for i, (z, x, y) in enumerate(all_points):
        component_dict[labels[i]].append((int(z), x, y))

    return component_dict

from sklearn.linear_model import LinearRegression

def fit_component_paths(component_dict, debug=False):
    """
    Fits linear functions x(z) and y(z) for each component.
    Returns a dict of {id: (model_x, model_y)}.
    If debug=True, plots scatter + fit lines.
    """
    fits = {}
    for cid, points in component_dict.items():
        zs = np.array([p[0] for p in points]).reshape(-1, 1)
        xs = np.array([p[1] for p in points])
        ys = np.array([p[2] for p in points])

        model_x = LinearRegression().fit(zs, xs)
        model_y = LinearRegression().fit(zs, ys)
        fits[cid] = (model_x, model_y)

    return fits

def get_rod_bounds_from_mask(rod_mask, margin=10):
    """
    Determine z-min and z-max from the extent of the rod mask.
    Trims 'margin' slices from both ends to avoid boundary artifacts.
    """
    z_indices = np.where(np.any(rod_mask, axis=(1, 2)))[0]

    if len(z_indices) < 2 * margin:
        raise ValueError("Not enough slices in rod mask to safely trim margins.")

    z_min = int(z_indices[0]) + margin
    z_max = int(z_indices[-1]) - margin
    return z_min, z_max

def create_rod_mask(shape, component_fits, z_min, z_max, radius=5, smoothed_img=None):
    Z, Y, X = shape

    # Compute average intensities if image is given
    intensity_list = []
    if smoothed_img is not None:
        for cid, (model_x, model_y) in component_fits.items():
            values = []
            for z in range(z_min, z_max + 1):
                cx = model_x.predict([[z]])[0]
                cy = model_y.predict([[z]])[0]
                yy, xx = np.ogrid[:Y, :X]
                circle = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
                values.extend(smoothed_img[z][circle])
            mean_intensity = np.mean(values)
            intensity_list.append((cid, mean_intensity))

        # Sort by intensity and re-map labels
        sorted_cids = [cid for cid, _ in sorted(intensity_list, key=lambda x: x[1])]
    else:
        # Keep original order if no intensities
        sorted_cids = sorted(component_fits.keys())

    final_mask = np.zeros((Z, Y, X), dtype=np.uint8)

    for new_label, cid in enumerate(sorted_cids, start=1):
        model_x, model_y = component_fits[cid]
        for z in range(z_min, z_max + 1):
            cx = model_x.predict([[z]])[0]
            cy = model_y.predict([[z]])[0]
            yy, xx = np.ogrid[:Y, :X]
            circle = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
            final_mask[z][circle] = new_label

    return final_mask


def relabel_by_model(labeled_img, model_name):
    """
    Remaps label values in `labeled_img` (1, 2, 3...) to model-specific values.

    Example:
        If model == "Mindways Model 3 CT", remaps 1→111, 2→112, ..., 5→115

    Parameters:
    - labeled_img (np.ndarray): 3D array with labels
    - model_name (str): Model string

    Returns:
    - relabeled_img (np.ndarray): 3D array with new labels
    """
    model_label_map = model_labels()

    if model_name not in model_label_map:
        raise ValueError(f"Model '{model_name}' not recognized.")
    
    print(f'Using Model: {model_name} for labelling')

    new_labels = model_label_map[model_name]
    unique_labels = sorted(np.unique(labeled_img[labeled_img > 0]))

    if len(unique_labels) != len(new_labels):
        raise ValueError(
            f"Label count mismatch: image has {len(unique_labels)} labels, model '{model_name}' expects {len(new_labels)}."
        )

    relabeled_img = np.zeros_like(labeled_img)
    for old_label, new_label in zip(unique_labels, new_labels):
        relabeled_img[labeled_img == old_label] = new_label

    return relabeled_img


def segment_phantom_rods(image_path, radius=10, model_name='Mindways Model 3', output_path=None):
    img_original = sitk.ReadImage(image_path)
    img_ras = sitk.DICOMOrient(img_original, 'RAS')

    spacing = [1.0, 1.0, 1.0]
    original_spacing = img_ras.GetSpacing()
    original_size = img_ras.GetSize()
    new_size = [int(round(osz * ospc / nspc)) for osz, ospc, nspc in zip(original_size, original_spacing, spacing)]

    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputDirection(img_ras.GetDirection())
    resampler.SetOutputOrigin(img_ras.GetOrigin())
    resampler.SetInterpolator(sitk.sitkLinear)
    img_resampled = resampler.Execute(img_ras)

    arr = sitk.GetArrayFromImage(img_resampled)
    smoothed = gaussian(arr, sigma=1.0, preserve_range=True)

    binary = smoothed > -250
    binary_filled = binary_fill_holes(binary)
    opened = binary_opening(binary_filled, footprint=np.ones((3, 3, 3)))

    labeled = label(opened)
    props = sorted(regionprops(labeled), key=lambda x: x.area, reverse=True)

    if len(props) < 2:
        raise RuntimeError("Could not identify the phantom holder.")

    holder_mask = labeled == props[1].label
    holder_mask = binary_erosion(holder_mask, iterations=2)

    mask_output = hough_detection(arr, holder_mask)

    centroids = extract_slicewise_centroids(mask_output)
    n_components = estimate_num_components(model_name)
    full_slices = filter_full_slices(centroids, n_components)
    component_dict = cluster_centroids(full_slices, n_components)
    component_fits = fit_component_paths(component_dict)
    z_min, z_max = get_rod_bounds_from_mask(mask_output, margin=10)

    final_mask = create_rod_mask(mask_output.shape, component_fits, z_min, z_max, radius=radius, smoothed_img=smoothed)
    relabeled_mask = relabel_by_model(final_mask, model_name)

    mask_ras = sitk.GetImageFromArray(relabeled_mask)
    mask_ras.CopyInformation(img_resampled)

    resampler.SetSize(img_ras.GetSize())
    resampler.SetOutputSpacing(img_ras.GetSpacing())
    resampler.SetOutputDirection(img_ras.GetDirection())
    resampler.SetOutputOrigin(img_ras.GetOrigin())
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    mask_orig_spacing = resampler.Execute(mask_ras)

    direction = img_original.GetDirection()
    mask_final = sitk.DICOMOrient(mask_orig_spacing, sitk.DICOMOrientImageFilter().GetOrientationFromDirectionCosines(direction))
    mask_final.CopyInformation(img_original)

    if output_path is not None:
        if os.path.isdir(output_path):
            base_name = os.path.basename(image_path).replace(".nii.gz", "").replace(".nii", "")
            output_file = os.path.join(output_path, f"{base_name}_calib_rods.nii.gz")
        else:
            output_file = output_path

        sitk.WriteImage(mask_final, output_file)
        print(f"Saved relabeled calibration rods to: {output_file}")

    return mask_final


def main():
    """
    Command-line interface for segmenting calibration rods in CT phantom images.

    This function parses input arguments and runs the segmentation pipeline using
    the provided CT image. It detects and labels rods based on the selected phantom model.

    Usage:
        ogoSegmentCalibrationRods <input.nii.gz> --output <output.nii.gz or directory>

    Arguments:
        <input.nii.gz>               Path to the input NIfTI image.
        --output                     Path to the output file or directory. If a directory is given,
                                     the output file will be named automatically.
        --model                      Name of the calibration phantom model to use for labeling.
                                     Default: "Mindways Model 3 CT"
        --radius                     Approximate radius of rods in mm. Default: 7

    Returns:
        None. Writes the labeled output segmentation to file if --output is provided.

    Matthias Walle
    """

    parser = argparse.ArgumentParser(description="Segment calibration rods in CT phantom images.")
    parser.add_argument("input", type=str, help="Path to input image (NIfTI). Supports shell expansion.")
    parser.add_argument("--output", type=str, default=None, help="Output path or folder.")
    parser.add_argument("--model", type=str, default="Mindways Model 3 CT", help="Name of the phantom model.")
    parser.add_argument("--radius", type=int, default=7, help="Rod radius in mm.")
    args = parser.parse_args()

    mask = segment_phantom_rods(
        image_path=args.input,
        radius=args.radius,
        model_name=args.model,
        output_path=args.output
    )


if __name__ == '__main__':
    main()
