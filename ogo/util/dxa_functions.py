# /------------------------------------------------------------------------------+
# | 27-March-2025                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+
#
# Functions that are used to generate 2D DXA ROIs from clinical CT scans. 
#
# 
#
import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from skimage.draw import circle_perimeter, disk



def load_image_and_mask(image_path, mask_path):
    image = sitk.ReadImage(image_path)
    mask = sitk.ReadImage(mask_path)
    return image, mask


def find_mask_bounds(mask_sitk):
    mask_array = sitk.GetArrayFromImage(mask_sitk)
    non_zero_indices = mask_array.nonzero()
    min_x, max_x = non_zero_indices[0].min(), non_zero_indices[0].max()
    min_y, max_y = non_zero_indices[1].min(), non_zero_indices[1].max()
    min_z, max_z = non_zero_indices[2].min(), non_zero_indices[2].max()
    return (min_x, max_x), (min_y, max_y), (min_z, max_z)

def crop_image(image_sitk, mask_sitk, buffer=30):

    #Crops an image and its corresponding mask image with a buffer on the outside so that it is not exactly down to the mask.
    
    (min_x, max_x), (min_y, max_y), (min_z, max_z) = find_mask_bounds(mask_sitk)
    
    min_x = max(min_x - buffer, 0)
    max_x = min(max_x + buffer, mask_sitk.GetSize()[2] - 1)
    min_y = max(min_y - buffer, 0)
    max_y = min(max_y + buffer, mask_sitk.GetSize()[1] - 1)
    min_z = max(min_z - buffer, 0)
    max_z = min(max_z + buffer, mask_sitk.GetSize()[0] - 1)
    
    extract_size = [int(max_z - min_z + 1), int(max_y - min_y + 1), int(max_x - min_x + 1)]
    extract_index = [int(min_z), int(min_y), int(min_x)]
    
    extractor = sitk.RegionOfInterestImageFilter()
    extractor.SetSize(extract_size)
    extractor.SetIndex(extract_index)
    
    cropped_image = extractor.Execute(image_sitk)
    cropped_mask = extractor.Execute(mask_sitk)
    
    return cropped_image, cropped_mask

def get_center(img):
    """
    This function returns the physical center point of a 3d sitk image
    :param img: The sitk image we are trying to find the center of
    :return: The physical center point of the image
    """
    width, height, depth = img.GetSize()
    return img.TransformIndexToPhysicalPoint((int(np.ceil(width/2)),
                                              int(np.ceil(height/2)),
                                              int(np.ceil(depth/2))))

def rotation3d(image, theta_z, interpolator, show=False):
    """
    This function rotates an image across each of the x, y, z axes by theta_x, theta_y, and theta_z degrees
    respectively
    :param image: An sitk MRI image
    :param theta_x: The amount of degrees the user wants the image rotated around the x axis
    :param theta_y: The amount of degrees the user wants the image rotated around the y axis
    :param theta_z: The amount of degrees the user wants the image rotated around the z axis
    :param interpolator: Type of interpolation (sitk.sitkNearestNeighbor, sitk.sitkLinear, sitk.sitkBSpline)
    :param show: Boolean, whether or not the user wants to see the result of the rotation
    :return: The rotated image
    """
    theta_z = np.deg2rad(theta_z)
    euler_transform = sitk.Euler3DTransform()
    image_center = get_center(image)
    euler_transform.SetCenter(image_center)

    direction = image.GetDirection()
    axis_angle = (direction[2], direction[5], direction[8], theta_z)
    np_rot_mat = matrix_from_axis_angle(axis_angle)
    euler_transform.SetMatrix(np_rot_mat.flatten().tolist())
    resampled_image = resample(image, euler_transform, interpolator)
    if show:
        print(euler_transform.GetMatrix())
        slice_num = int(input("Enter the index of the slice you would like to see"))
        plt.imshow(sitk.GetArrayFromImage(resampled_image)[slice_num])
        plt.show()
    return resampled_image

def matrix_from_axis_angle(a):
    """ Compute rotation matrix from axis-angle.
    This is called exponential map or Rodrigues' formula.
    Parameters
    ----------
    a : array-like, shape (4,)
        Axis of rotation and rotation angle: (x, y, z, angle)
    Returns
    -------
    R : array-like, shape (3, 3)
        Rotation matrix
    """
    ux, uy, uz, theta = a
    c = np.cos(theta)
    s = np.sin(theta)
    ci = 1.0 - c
    R = np.array([[ci * ux * ux + c,
                   ci * ux * uy - uz * s,
                   ci * ux * uz + uy * s],
                  [ci * uy * ux + uz * s,
                   ci * uy * uy + c,
                   ci * uy * uz - ux * s],
                  [ci * uz * ux - uy * s,
                   ci * uz * uy + ux * s,
                   ci * uz * uz + c],
                  ])

    # This is equivalent to
    # R = (np.eye(3) * np.cos(theta) +
    #      (1.0 - np.cos(theta)) * a[:3, np.newaxis].dot(a[np.newaxis, :3]) +
    #      cross_product_matrix(a[:3]) * np.sin(theta))

    return R


def transform_mask(mask, rotation, label):
    mask_transformed = rotation3d(mask, rotation, sitk.sitkNearestNeighbor, False)        
    maskProjectionFilt = sitk.BinaryProjectionImageFilter() 
    seg_2d = maskProjectionFilt.Execute(mask)
    seg_2d = seg_2d<0
    mask = mask_transformed==label
    mask_projection = maskProjectionFilt.Execute(mask)
    mask_projection = label * mask_projection
    bin_seg_2d = seg_2d>0
    bin_part_2d = mask_projection>0
    overlap_2d = (bin_seg_2d + bin_part_2d)==2
    mask_2d = 1 - overlap_2d
    this_label_2d = sitk.Mask(mask_projection, mask_2d)
    seg_2d = seg_2d + label*(this_label_2d>0)
    return seg_2d

# Adapted from: https://www.geeksforgeeks.org/bresenhams-algorithm-for-3-d-line-drawing/ 
def bresenham_line_3d(start, end):
    points = []

    # Calculate differences between start and stop 
    dx = abs(end[0] - start[0])
    dy = abs(end[1] - start[1])
    dz = abs(end[2] - start[2])
    
    # calculating direction of steps
    xs = 1 if end[0] > start[0] else -1
    ys = 1 if end[1] > start[1] else -1
    zs = 1 if end[2] > start[2] else -1

    # Driving axis is X-axis
    if dx >= dy and dx >= dz:
        p1 = 2 * dy - dx
        p2 = 2 * dz - dx
        while start[0] != end[0]:
            start[0] += xs
            if p1 >= 0:
                start[1] += ys
                p1 -= 2 * dx
            if p2 >= 0:
                start[2] += zs
                p2 -= 2 * dx
            p1 += 2 * dy
            p2 += 2 * dz
            points.append((start[0], start[1], start[2]))
            
    # Driving axis is Y-axis
    elif dy >= dx and dy >= dz:
        p1 = 2 * dx - dy
        p2 = 2 * dz - dy
        while start[1] != end[1]:
            start[1] += ys
            if p1 >= 0:
                start[0] += xs
                p1 -= 2 * dy
            if p2 >= 0:
                start[2] += zs
                p2 -= 2 * dy
            p1 += 2 * dx
            p2 += 2 * dz
            points.append((start[0], start[1], start[2]))
            
   # Driving axis is Z-axis
    else:
        p1 = 2 * dy - dz
        p2 = 2 * dx - dz
        while start[2] != end[2]:
            start[2] += zs
            if p1 >= 0:
                start[1] += ys
                p1 -= 2 * dz
            if p2 >= 0:
                start[0] += xs
                p2 -= 2 * dz
            p1 += 2 * dy
            p2 += 2 * dx
            points.append((start[0], start[1], start[2]))

    return points

# Adapted from: https://www.geeksforgeeks.org/bresenhams-algorithm-for-3-d-line-drawing/ 
def bresenham_line_2d(start, end):
    points = []

    # Calculate differences between start and stop 
    dx = abs(end[0] - start[0])
    dy = abs(end[1] - start[1])
    
    # calculating direction of steps
    xs = 1 if end[0] > start[0] else -1
    ys = 1 if end[1] > start[1] else -1

    # Driving axis is X-axis
    if dx > dy:
        p = 2 * dy - dx
        while start[0] != end[0]:
            start[0] += xs
            if p >= 0:
                start[1] += ys
                p -= 2 * dx
            p += 2 * dy
            points.append((start[0], start[1]))
    # Driving axis is Y-axis
    else:
        p = 2 * dx - dy
        while start[1] != end[1]:
            start[1] += ys
            if p >= 0:
                start[0] += xs
                p -= 2 * dy
            p += 2 * dx
            points.append((start[0], start[1]))

    return points

def plot_hough_transform(cx,cy,radii,img):
    mask = np.zeros_like(img, dtype=float)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    for center_y, center_x, radius in zip(cy, cx, radii):
        circy, circx = disk((center_y, center_x), radius, shape=img.shape)
        mask[circy, circx] = 1

    ax.imshow(img, cmap=plt.cm.gray)
    ax.imshow(mask, cmap="Reds", alpha=mask/2)
    plt.show()
    

def calculate_angle_between_lines(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    angle_rad = np.arccos(dot_product / (norm_v1 * norm_v2))
    angle_deg = np.degrees(angle_rad)
    return angle_deg