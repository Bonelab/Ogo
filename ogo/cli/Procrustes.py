# /------------------------------------------------------------------------------+
# | 24-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import numpy as np
import SimpleITK as sitk
from scipy.spatial import procrustes
import matplotlib.pyplot as plt

from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb


def Procrustes(input_filenames, label_list, sensitivity, plot):
    # Master list of labels
    labelsDict = lb.labels_dict

    # Check input files exist and correct format
    ogo.message('Reading input files:')
    for input_image in input_filenames:
        if not os.path.isfile(input_image):
            os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))
        if not (input_image.lower().endswith('.nii') or input_image.lower().endswith('.nii.gz')):
            os.sys.exit('[ERROR] Input must be type NIFTI file: \"{}\"'.format(input_image))
        ogo.message('  {}'.format(input_image))

    # Check minimum number of input images
    if len(input_filenames) < 2:
        os.sys.exit('[ERROR] At least two files are required for Procrustes analysis.')

    base = True  # Used to identify the base case for Procrustes
    disparity = 0.0

    # Stuff for the plot
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Label Positions")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")

    # Calculate centroids for each label
    for input_image in input_filenames:
        fname = os.path.basename(input_image)
        arr = []

        ct = sitk.ReadImage(input_image)
        filt = sitk.LabelShapeStatisticsImageFilter()
        filt.Execute(ct)
        n_labels = filt.GetNumberOfLabels()
        labels = filt.GetLabels()

        if label_list is None:
            label_list = labels

        print()
        print('{}'.format(fname))
        print('  {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:8s}'.format('Label', 'Voxels', 'Cx', 'Cy', 'Cz', 'Name'))
        for idx in range(n_labels):
            label = labels[idx]
            if label in label_list:
                centroid = filt.GetCentroid(label)
                arr.append(centroid)
                n_pixels = filt.GetNumberOfPixels(label)
                print('  {:8d} {:8d} {:8.3f} {:8.3f} {:8.3f} {}'.format(label, n_pixels, centroid[0], centroid[1],
                                                                        centroid[2], labelsDict[label]['LABEL']))

        # Procrustes
        if base:
            arr_base = np.array(arr)
            base = False
            # scat_plot = ax.scatter(mtx1[:,0],mtx1[:,1],mtx1[:,2])

        else:
            arr_curr = np.array(arr)
            mtx1, mtx2, disparity = procrustes(arr_base, arr_curr)
            print('   DISPARITY = {:8.6f}'.format(disparity))
            if (disparity > sensitivity):
                print('   WARNING!!! There is a possible blunder in the voxel labels for this dataset.')
            scat_plot = ax.scatter(mtx2[:, 0], mtx2[:, 1], mtx2[:, 2])
            scat_plot = ax.scatter(mtx1[:, 0], mtx1[:, 1], mtx1[:, 2], marker="x", c="red")

    if plot:
        plt.show()

    # Plotting
    # https://likegeeks.com/3d-plotting-in-python/
    # c = np.array([[1, 3, 2], [1, 2, 1], [1, 1, 7], [2, 1, 0]], 'd')
    # fig = plt.figure(figsize=(4,4))
    # ax = fig.add_subplot(111, projection='3d')


#    ax.scatter(2,3,4) # plot the point (2,3,4) on the figure

# np.random.seed(42)
# ages = np.random.randint(low = 8, high = 30, size=35)
# xs = arr_base[:,0]
# ys = arr_base[:,1]
# zs = arr_base[:,2]
# weights = np.random.randint(30, 160, 35)
# gender_labels = np.random.choice([0, 1], 35) #0 for male, 1 for female


#    ax.scatter(xs,ys,zs, marker="x", c="red")
#    ax.set_xlim(100,200)
#    ax.set_ylim(20,160)
#    ax.set_zlim(5,35)
#    ax.set_xticks([100,125,150,175,200])

# ax.scatter(c[0,:],c[1,:],c[2,:])
# scat_plot = ax.scatter(arr_base[:,0],arr_base[:,1],arr_base[:,2])
# scat_plot = ax.scatter(arr_curr[:,0],arr_curr[:,1],arr_curr[:,2])
# scat_plot = ax.scatter(xs = heights, ys = weights, zs = ages, c=gender_labels)
# cb = plt.colorbar(scat_plot, pad=0.2)
# cb.set_ticks([0,1])
# cb.set_ticklabels(["Male", "Female"])
# ax.scatter(xs = heights, ys = weights, zs = ages, c=gender_labels)

# labels = ["Florida", "Georgia", "California"]
#
# for l in labels:
#    ages = np.random.randint(low = 8, high = 20, size=20)
#    heights = np.random.randint(130, 195, 20)
#    weights = np.random.randint(30, 160, 20)
#    ax.scatter(xs = heights, ys = weights, zs = ages, label=l)
#
# ax.set_title("Age-wise body weight-height distribution")
# ax.set_xlabel("Height (cm)")
# ax.set_ylabel("Weight (kg)")
# ax.set_zlabel("Age (years)")
# ax.legend(loc="best")


def main():
    # Setup description
    description = '''
Calculates the disparity between labelled datasets using Procrustes
algorithm. This is usesful for finding bones that are incorrectly
labelled. 

If an L1 and L2 are swapped, or L4 and L5 are swapped, or left and
right femur are swapped then a sensitivity of 0.01 is sufficient.

'''
    epilog = '''
Example calls: 
ogoProcrustes mask1.nii mask2.nii mask3.nii
ogoProcrustes mask1.nii mask2.nii mask3.nii --plot
ogoProcrustes mask1.nii mask2.nii mask3.nii --label_list 1 2 3 4 5 6 7 8 9 10
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoProcrustes",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_filenames', nargs='*', help='Input mask image files (*.nii, *.nii.gz)')
    parser.add_argument('--label_list', type=int, nargs='*', default=None, metavar='ID',
                        help='Specify which labels to analyze; space separated (e.g. 1 2 3) (default: all)')
    parser.add_argument('--sensitivity', type=float, default=0.01, metavar='VAL',
                        help='Sensitivity to trigger warning (default: %(default)s)')
    parser.add_argument('--plot', action='store_true', help='Plot the label points')

    print()

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('Procrustes', vars(args)))

    # Run program
    Procrustes(**vars(args))


if __name__ == '__main__':
    main()
