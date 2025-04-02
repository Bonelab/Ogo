import sys
import os
import argparse
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy.ndimage import label
import scipy.ndimage as ndimage
from sklearn.cluster import DBSCAN
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.draw import circle_perimeter, disk
from scipy.ndimage.measurements import label
import matplotlib.pyplot as plt
from skimage.draw import polygon
from scipy.ndimage import binary_fill_holes


import ogo.util.dxa_functions as ogo
import ogo.util.Helper as helper
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_txt import write_txt


def GenerateDXAROI(image_filename, mask_filename, output_path_image, output_path_mask, region, aBMD, overwrite=False):
    helper.message("Starting GenerateDXAROI script...")
    # Check if output exists and should overwrite
    if os.path.isfile(output_path_mask) and not overwrite:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_path_image))
        if result.lower() not in ['y', 'yes']:
            helper.message('Not overwriting. Exiting...')
            os.sys.exit()
    
    # Read input image and mask 
    ct, mask = ogo.load_image_and_mask(image_filename, mask_filename)
    spacing_unprojected = ct.GetSpacing()

    #cropping ct and mask down to only the bone of interest 
    cropped_ct, cropped_mask = ogo.crop_image(ct, mask)

    #projecting the image and mask down into 2D 
    projectionFilt = sitk.SumProjectionImageFilter()
    projectionFilt.SetProjectionDimension(1)
    ct_projection = projectionFilt.Execute(cropped_ct)
    sitk.WriteImage(ct_projection, output_path_image)

    mask_projectionFilt = sitk.BinaryProjectionImageFilter() 
    mask_projectionFilt.SetProjectionDimension(1)
    mask_projection = mask_projectionFilt.Execute(cropped_mask)
    sitk.WriteImage(mask_projection, output_path_mask)

    #masking the projection ... 
    masked_projection = sitk.Mask(ct_projection, mask_projection)
    spacing_projected = masked_projection.GetSpacing()

    #converting the numpy for analysis 
    ct_projection_np = sitk.GetArrayFromImage(masked_projection)
    mask_projection_np = sitk.GetArrayFromImage(mask_projection)
    masked_projection_np = sitk.GetArrayFromImage(masked_projection)

    ct_numpy = sitk.GetArrayFromImage(ct)
    mask_numpy = sitk.GetArrayFromImage(mask)
    
    if region == "FE":
        ct_numpy = sitk.GetArrayFromImage(ct)
        mask_numpy = sitk.GetArrayFromImage(mask)  
        highest_indices = np.argmax(mask_numpy, axis=0)
        topmost_indices = np.where(mask_numpy[highest_indices, 
                                              np.arange(mask_numpy.shape[1])[:, None], np.arange(mask_numpy.shape[2])], 1, 0)
        cutoff_index = np.max(highest_indices) - 120 #length of cut off  

        for x in range(mask_numpy.shape[1]):  
            for z in range(mask_numpy.shape[2]):  
                mask_numpy[cutoff_index > np.arange(mask_numpy.shape[0]), x, z] = 0  
        helper.message('Writing out mask for FE to: {}'.format(output_path_mask))
        helper.message('WARNING: this mask will only work on the original full sized image (not the cropped image...)')
        modified_sitk_mask = sitk.GetImageFromArray(mask_numpy)
        sitk.WriteImage(modified_sitk_mask, output_path_mask) # NOTE: this is for the full image (not the cropped one...)
    

    #cropping 2D projection... 
    last_row_with_1 = -1 
    for i in range(mask_projection_np.shape[0] - 1, -1, -1):
        if 1 in mask_projection_np[i]:
            last_row_with_1 = i
            break
    cutoff_length = 130
    end_row = last_row_with_1 - cutoff_length
    for i in range(end_row):
        for j in range(mask_projection_np.shape[2]):
            if mask_projection_np[i,0,j] == 1:
                mask_projection_np[i,0,j] = 0

    shaft_array = []
    lengths = []
    for i in range(end_row, end_row+5, 1):
        for j in mask_projection_np[i,0,:]:
            if j > 0: 
                shaft_array.append(j)
        length = len(shaft_array)
        shaft_array = []
        lengths.append(length)
    lengths = np.array(lengths)

    halves = []
    for length in lengths:
        half = math.ceil(length / 2)
        halves.append(half)
    halves = np.array(halves)
    mean_column = np.mean(halves)

    # finding the centerline of the femoral shaft
    for idx, value in enumerate(mask_projection_np[end_row,0,:]):
        if value > 0:
            new_index = int(idx + mean_column)
            if new_index < mask_projection_np.shape[2]:  
                for row_idx in range(end_row, mask_projection_np.shape[0],1):
                    if mask_projection_np[row_idx, 0, new_index,] > 0:
                        mask_projection_np[row_idx, 0,new_index] = 2
                    else:
                        break 
            break  
    
    # identifying the medial and lateral sides of the femur
    for row_idx in range(end_row, mask_projection_np.shape[0],1):
        first_index = -1
        last_index = -1
        for idx, value in enumerate(mask_projection_np[row_idx, 0, :]):
            if value > 0:
                if first_index == -1:
                    first_index = idx
                last_index = idx
        # Change the values at the first and last index to 3
        if first_index != -1:
            mask_projection_np[row_idx, 0, first_index] = 3
        if last_index != -1 and last_index != first_index:
            mask_projection_np[row_idx, 0,last_index ] = 4

    
    # looking through the medial side of the femur to find the lesser trochanter
    edge_pixels = np.argwhere(mask_projection_np == 3)
    edge_columns = edge_pixels[:, 2] 
    edge_rows = edge_pixels[:, 0] 
    locations = np.vstack((edge_rows, edge_columns)).T
    intensities = []

    for i,j in locations: 
        intensity = ct_projection_np[i,0,j]
        intensities.append(intensity)
    
    #location in image where mask == 3 in form of [row column intensity]
    locations_with_intensities = np.vstack((edge_rows, edge_columns, intensities)).T

    #taking the spatial derivative of the column (higher derivative == sharp change in column location)
    spatial_derivative = np.gradient(locations_with_intensities[:,1])
    negative_derivative_indices = np.argwhere(spatial_derivative < 0).flatten()

    #changes all pixels with negative spatial derivative to 5 (identifying potential lesser trochanter locations)
    original_indices = edge_pixels[negative_derivative_indices]
    for i in original_indices[:len(original_indices)//2]:
        mask_projection_np[i[0],0, i[2]] = 5

    #This code identifies all clusters of pixels with a value of 5 in mask_projection_np using DBSCAN.
    #It then finds the largest cluster and changes the pixels in that cluster to 6 while setting 
    # the rest of the pixels with value 5 to 0.
    coords = np.column_stack(np.where(mask_projection_np == 5))
    db = DBSCAN(eps=5, min_samples=1).fit(coords)
    labels = db.labels_
    unique_labels, counts = np.unique(labels, return_counts=True)
    largest_cluster_label = unique_labels[np.argmax(counts)]
    temp_mask = np.zeros_like(mask_projection_np)
    largest_cluster_coords = coords[labels == largest_cluster_label]
    temp_mask[largest_cluster_coords[:, 0], 0, largest_cluster_coords[:, 2]] = 6
    mask_projection_np[temp_mask == 6] = 5

    #finding spatial location of the lesser trochanter
    trochanter_spots = []
    for i in range(mask_projection_np.shape[0]):
        for j in range(mask_projection_np.shape[2] - 1):
            if mask_projection_np[i,0,:][j] == 5:
                trochanter_spots.append(i)
                
    trochanter_start = trochanter_spots[0]

    #turning anything above where the trochanter starts to 0 to cut off the mask below the lesser trochanters
    for i in range(trochanter_start):
        for j in range(mask_projection_np.shape[2] - 1):
            if mask_projection_np[i,0,:][j] > 0:
                mask_projection_np[i,0,:][j] = 0
    
    #turning all the different regions to value 1 (binary mask)
    mask_projection_np[mask_projection_np == 2] = 1
    mask_projection_np[mask_projection_np == 6] = 1
    mask_projection_np[mask_projection_np == 3] = 1
    mask_projection_np[mask_projection_np == 4] = 1
    mask_projection_np[mask_projection_np == 5] = 1

    total_hip_mask = sitk.GetImageFromArray(mask_projection_np)
    masked_projection_np[:,0,:][mask_projection_np[:,0,:] != 1] = 0
    sitk_masked_projection = sitk.GetImageFromArray(masked_projection_np)

    #Canny edge detection (for Hough transform)
    edges = sitk.CannyEdgeDetection(sitk.Cast(sitk_masked_projection, sitk.sitkFloat32), 
                                lowerThreshold=0.0, upperThreshold=100.0, variance = (5.0,5.0,5.0))
    mask_of_edges = sitk.CannyEdgeDetection(sitk.Cast(total_hip_mask, sitk.sitkFloat32), 
                                lowerThreshold=0.0, upperThreshold=0.7)
    mask_of_edges = sitk.BinaryDilate(sitk.Cast(mask_of_edges, sitk.sitkUInt8), [2,2,2])
    edges_casted = sitk.Cast(edges, sitk.sitkUInt8)
    mask_of_edges_casted = sitk.Cast(mask_of_edges, sitk.sitkUInt8)
    #creating a mask that is just the outer edges of the total hip mask (don't care about edges within the mask)
    combined_mask = sitk.And(edges_casted, mask_of_edges_casted)

    #Hough transform to find the femoral head
    combined_mask_np = sitk.GetArrayFromImage(combined_mask)
    hough_radius_to_try = np.arange(32,40,30)
    hough_res = hough_circle(combined_mask_np[:,0,:], hough_radius_to_try)
    accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radius_to_try, total_num_peaks=1)

    total_hip_projected_np = sitk.GetArrayFromImage(sitk_masked_projection)
    #img = total_hip_projected_np[:,0,:]

    y_center = cx[0]
    x_center = cy[0]
    rad = radii[0]

    # for the left femur (potentially not the same for the right femur)
    #In hologic scanners, the top of the femoral neck begins at the two points on the left and right of the femoral neck 
    #that are the smallest distance to the center of the femoral head 
    # so the following code is finding the distance of each point on the edge to the center of the femoral neck and 
    # storing the distances so that it can keep the two shortest points. 
    # it stores their coordinates in coordinates_left and coordinates_right. 
    distances_left = []
    positions_left = []
    distances_right = []
    positions_right = []

    combined_mask_np = combined_mask_np[:,0,:]

    for x in range(combined_mask_np.shape[0]):
        for y in range(combined_mask_np.shape[1]):
            if combined_mask_np[x, y] == 1 and x < x_center and y < y_center:
                distance_left = np.sqrt((x - x_center)**2 + (y - y_center)**2)
                distances_left.append(distance_left)
                positions_left.append((x,y))
                
    minimum_distance_left = np.min(distances_left)
    minimum_distance_index_left = np.where(distances_left == minimum_distance_left)
    coordinates_left = positions_left[minimum_distance_index_left[0][0]]
        
    adjusted = x_center + (rad/2)
    for x in range(combined_mask_np.shape[0]):
        for y in range(combined_mask_np.shape[1]):
            if combined_mask_np[x, y] == 1 and y > y_center and x > coordinates_left[0] and x < adjusted:
                distance_right = np.sqrt((x - x_center)**2 + (y - y_center)**2)
                distances_right.append(distance_right)
                positions_right.append((x,y))
                
    minimum_distance_right = np.min(distances_right)
    minimum_distance_index_right = np.where(distances_right == minimum_distance_right)
    coordinates_right = positions_right[minimum_distance_index_right[0][0]]

    edges_np = sitk.GetArrayFromImage(edges)
    connected_lines = np.zeros_like(edges_np[:,0,:])
    # turning the closest points to 1
    connected_lines[coordinates_left[0],coordinates_left[1]] = 1
    connected_lines[coordinates_right[0],coordinates_right[1]] = 1
    
    
     #function for the left femur, creating the "edges" of the femoral neck 
    for j in range(connected_lines.shape[0]):
        for i in range(connected_lines.shape[1]):
            if connected_lines[j, i] == 1 and combined_mask_np[j, i] == 1:
                for dj in range(-10,0):
                    for di in range(0, 10):
                        new_j, new_i = j + dj, i + di
                        if 0 <= new_j < connected_lines.shape[0] and 0 <= new_i < connected_lines.shape[1]:
                            if combined_mask_np[new_j, new_i] == 1:
                                connected_lines[new_j, new_i] = 1

    #relabelling the left and right sides to labels 1 and 2
    structure = np.ones((3, 3))
    labeled, ncomponents = label(connected_lines, structure)
    line1_points = np.argwhere(labeled == 1)
    line2_points = np.argwhere(labeled == 2)

    #finding the narrowest point of the femoral neck (the two points on the left and right that are closest together...)
    min_length = min(len(line1_points), len(line2_points))

    distances = []
    closest_index_line1 = None
    closest_index_line2 = None
    min_distance = float('inf')

    for i in range(min_length):
        for j in range(min_length):
            point1 = line1_points[i]
            point2 = line2_points[j]

            distance = np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
            distances.append(distance)

            if distance < min_distance:
                min_distance = distance
                closest_index_line1 = i
                closest_index_line2 = j

    labeled2 = np.zeros_like(labeled)
    if closest_index_line1 is not None:
        labeled2[tuple(line1_points[closest_index_line1])] = 1
        labeled2[tuple(line2_points[closest_index_line2])] = 1

    center = np.array([x_center, y_center])
    x1 = tuple(line1_points[closest_index_line1])[0]
    x2 = tuple(line2_points[closest_index_line2])[0]
    y1 = tuple(line1_points[closest_index_line1])[1]
    y2 = tuple(line2_points[closest_index_line2])[1]

    # connecting the two narrowest points (line)
    line_points_narrowest = ogo.bresenham_line_2d(line1_points[closest_index_line1], line2_points[closest_index_line2])
    femoral_neck_narrowest = np.zeros_like(labeled)
    for point in line_points_narrowest:
        femoral_neck_narrowest[tuple(point)] = 1
        
    femoral_neck_narrowest[x_center, y_center] = 1

    #finding the point on the narrowest line that creates an angle closest to 90 degrees with the center of the femoral head. 
    #This intersection of the line + center of the femoral head becomes the axis for the femoral neck ROI aka the (0,0)
    best_point = None
    closest_angle = float('inf')

    for points in line_points_narrowest:
        line2_vector_end = np.array([points[0], points[1]])
        vector_line2 = line2_vector_end - center
        
        vector_line1 = np.array([x2 - x1, y2 - y1])
        
        
        angle = ogo.calculate_angle_between_lines(vector_line1, vector_line2)
        
        if abs(angle - 90) < abs(closest_angle - 90):
            closest_angle = angle
            best_point = points

    helper.message('The best point on the narrowest point on the femoral neck is: {}'.format(best_point))
    helper.message('The closest angle to 90 degrees is: {}'.format(closest_angle))

    line_point_midpoint = ogo.bresenham_line_2d(center, best_point)
    for points in line_point_midpoint:
        femoral_neck_narrowest[tuple(points)] = 1

    #creating femoral neck box
    best_point_x, best_point_y = best_point[0], best_point[1]
    width, height = 82, 34

    x_axis_dir = np.array([x2 - x1, y2 - y1])
    x_axis_dir = x_axis_dir / np.linalg.norm(x_axis_dir)

    y_axis_dir = np.array([best_point_x - x_center, best_point_y - y_center])
    y_axis_dir = y_axis_dir / np.linalg.norm(y_axis_dir)

    half_width = width // 2
    half_height = height // 2
    fourth_height = height // 4
    sixth_height = height // 6

    top_left = np.array([best_point_x, best_point_y]) - half_width * x_axis_dir + sixth_height * x_axis_dir
    top_right = np.array([best_point_x, best_point_y]) + half_width * x_axis_dir + sixth_height * x_axis_dir
    bottom_left = np.array([best_point_x, best_point_y]) - half_width * x_axis_dir - half_height * y_axis_dir
    bottom_right = np.array([best_point_x, best_point_y]) + half_width * x_axis_dir - half_height * y_axis_dir

    corners = np.array([top_left, top_right, bottom_right, bottom_left], dtype=int)

    femoral_neck_points = np.zeros_like(connected_lines)
    femoral_neck_box = np.zeros_like(femoral_neck_narrowest)

    # filling in the femoral neck (with 1's)
    rr, cc = polygon(corners[:, 0], corners[:, 1], femoral_neck_box.shape)
    femoral_neck_box[rr, cc] = 1
    mask_np = sitk.GetArrayFromImage(mask_projection)
    added = np.add(mask_np[:,0,:], femoral_neck_box)
    femoral_neck_region = np.zeros_like(femoral_neck_box)
    femoral_neck_region[added == 2] = 1
    femoral_neck_mask = np.zeros_like(masked_projection_np)
    femoral_neck_mask[:,0,:] = femoral_neck_region
    femoral_neck_mask = sitk.GetImageFromArray(femoral_neck_mask)
    femoral_neck_mask.CopyInformation(masked_projection)
    if region == 'femoral_neck':
        helper.message('Writing out femoral neck mask to: {}'.format(output_path_mask))
        region = femoral_neck_mask
        sitk.WriteImage(femoral_neck_mask, output_path_mask)

    # isolating the top of the femoral neck, to subtract from the entire hip segmentation to just get the
    # total hip
    line_point_top = ogo.bresenham_line_2d(corners[2], corners[3])
    for points in line_point_top:
        femoral_neck_points[tuple(points)] = 1

    mask_without_line = mask_projection_np[:,0,:] * (1 - femoral_neck_points) #avoids having negatives
    labeled_mask, num_features = ndimage.label(mask_without_line)
    component_sizes = ndimage.sum(mask_without_line, labeled_mask, range(1, num_features + 1))

    sorted_indices = np.argsort(component_sizes)[::-1]  
    largest_component_label = sorted_indices[0] + 1 
    second_largest_component_label = sorted_indices[1] + 1 

    total_hip_region = (labeled_mask == largest_component_label).astype(int)
    femoral_head_region = (labeled_mask == second_largest_component_label).astype(int)

    total_hip_mask = np.zeros_like(mask_np)
    total_hip_mask[:,0,:] = total_hip_region
    if region == "total_hip":
        helper.message('Writing out total hip mask to: {}'.format(output_path_mask))
        total_hip_mask = sitk.GetImageFromArray(total_hip_mask)
        region = total_hip_mask
        sitk.WriteImage(total_hip_mask, output_path_mask)
    
    if aBMD: 
        region_masked = sitk.Mask(ct_projection, sitk.Cast(region, sitk.sitkUInt8))
        region_masked_np = sitk.GetArrayFromImage(region_masked)

        slices = cropped_ct.GetSize()
        y_slices = slices[1]

        mean_projected_vBMD_image = region_masked_np / y_slices

        #calculating voxel volume
        voxel_volume_mm3 = spacing_projected[0] * spacing_projected[1] * spacing_projected[2]
        voxel_volume_cm3 = voxel_volume_mm3 / 1000

        #calculating pixel area
        pixel_area_mm2 = spacing_projected[0] * spacing_projected[2]
        pixel_area_cm2 = pixel_area_mm2 / 100

        #turning units from mg/ccm to mg  
        image_mg_k2hpo4 = mean_projected_vBMD_image[:,0,:] * voxel_volume_cm3
        image_g_k2hpo4 = image_mg_k2hpo4 * 0.001

        #calculating total area under mask (cm^2)
        binary_mask = femoral_neck_region > 0
        total_area = binary_mask.sum() * pixel_area_cm2

        #calculating BMC
        BMC = image_g_k2hpo4.sum()

        #calculating aBMD (g/cm^2)
        aBMD_avg = BMC / total_area

        helper.message('Total area: {:8.4g}'.format(total_area))
        helper.message('Bone mineral content (BMC): {:8.4g}'.format(BMC))
        helper.message('Bone mineral density (BMD): {:8.4g}'.format(aBMD_avg))

    sys.exit()
    
    
    


def main():
    description = '''
    This script will calculate DXA regions of interest (ROIs) from clincial CT scans. Currently this is only in 2D: it takes in a 
    3D clinical CT scan and will output a 2D projection of the region of interest (i.e., hip, spine) and a mask of its related
    DXA ROI. 


    The DXA ROIs in the hip are: 

    1. Total hip 
    2. Femoral neck 

    And in the spine it simply projects downwards, to mimic the DXA scan. 

    WARNING: If your input image is not calibrated then the results here will be 
    incorrect. There is no calibration done as part of this application.

    EVEN BIGGER WARNING: This is currently still in the testing phase so use with caution! 

'''
    epilog = '''
    Example calls:

    ogoGenerateDXAROI image_k2hpo4.nii mask.nii.gz

    '''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoGenerateDXAROI",
        description=description,
        epilog=epilog
    )

    parser.add_argument('image_filename', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('mask_filename', help='Input image mask file (*.nii, *.nii.gz)')
    parser.add_argument('output_path_image',
                        help='Path to where the output 2D projected image will be output')
    parser.add_argument('output_path_mask', help='Output file name for 2D projected ROI mask')
    parser.add_argument('--region', default='femoral_neck',
                                choices=['femoral_neck', 'total_hip', 'L1', 'L2',
                                         'L3', 'L4', 'L1-L4', 'FE'],
                                help='Specify which DXA ROI (default: %(default)s)')
    parser.add_argument('--aBMD', action='store_true', 
                        help='Use when you would also like to evaluate the aBMD under the mask.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')
    print()

    # Parse and display
    args = parser.parse_args()

    if args.output_path_image:
        print(echo_arguments('GenerateDXAROI', vars(args)))

    # Run program
    GenerateDXAROI(**vars(args))


if __name__ == '__main__':
    main()
