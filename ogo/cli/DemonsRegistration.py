# /------------------------------------------------------------------------------+
# | 20-JUN-2023                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+ 


#imports
import SimpleITK as sitk 
import numpy as np
import argparse
import csv 
import matplotlib.pyplot as plt
from ogo.util.echo_arguments import echo_arguments


def DemonsRegistration(atlas_path, atlas_mask_path, new_image, output_image, output_graph, iterations, update_field, smooth_displacement):
    # read in the image
    moving_image = sitk.ReadImage(atlas_path)
    moving_mask = sitk.ReadImage(atlas_mask_path)
    fixed_image = sitk.ReadImage(new_image)

    # Gaussian filter the images
    moving_image = sitk.SmoothingRecursiveGaussian(moving_image)
    fixed_image = sitk.SmoothingRecursiveGaussian(fixed_image)


    def smoothing_and_resampling(image, shrink_factors, smoothing_sigmas):
    # adapted from http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/66_Registration_Demons.html
        if np.isscalar(shrink_factors):
            shrink_factors = [shrink_factors]*image.GetDimension()
        if np.isscalar(smoothing_sigmas):
            smoothing_sigmas = [smoothing_sigmas]*image.GetDimension()

        smoothed_image = sitk.SmoothingRecursiveGaussian(image, smoothing_sigmas)
        
        original_spacing = image.GetSpacing()
        original_size = image.GetSize()
        new_size = [int(sz/float(sf) + 0.5) for sf,sz in zip(shrink_factors,original_size)]
        new_spacing = [((original_sz-1)*original_spc)/(new_sz-1) 
                    for original_sz, original_spc, new_sz in zip(original_size, original_spacing, new_size)]
        return sitk.Resample(smoothed_image, new_size, sitk.Transform(), 
                            sitk.sitkLinear, image.GetOrigin(),
                            new_spacing, image.GetDirection(), 0.0, 
                            image.GetPixelID())

    def multiscale_demons(registration_algorithm,
                        fixed_image, moving_image, initial_transform = None, 
                        shrink_factors=None, smoothing_sigmas=None):

        if initial_transform:
            initial_displacement_field = sitk.TransformToDisplacementField(
                initial_transform, sitk.sitkVectorFloat64, fixed_image.GetSize(),
                fixed_image.GetOrigin(), fixed_image.GetSpacing(), fixed_image.GetDirection())

        else:
            initial_displacement_field = sitk.Image(fixed_image.GetWidth(), 
                                                    fixed_image.GetHeight(),
                                                    fixed_image.GetDepth(),
                                                    sitk.sitkVectorFloat64)
            initial_displacement_field.CopyInformation(fixed_image)

        # the multiscale part of the registration 
        if shrink_factors:
            for i in range(len(shrink_factors)):
                    resampled_fixed_image = smoothing_and_resampling(fixed_image, shrink_factors[i], smoothing_sigmas[i])
                    resampled_moving_image = smoothing_and_resampling(moving_image, shrink_factors[i], smoothing_sigmas[i])
                    initial_displacement_field = sitk.Resample (initial_displacement_field, resampled_fixed_image)
                    initial_displacement_field = registration_algorithm.Execute(resampled_fixed_image, resampled_moving_image, initial_displacement_field)
        
        
        full_res =  registration_algorithm.Execute(fixed_image, moving_image, sitk.Resample(initial_displacement_field, fixed_image))
        
        return full_res


    # creating callback to track registration metrics 
    metrics = []
    def iteration_callback(filter):
        print('\r{0}: {1:.2f}'.format(filter.GetElapsedIterations(), filter.GetMetric()), end='')
        metrics.append(filter.GetMetric())

    # setting the initial transform
    initial_transform = sitk.CenteredTransformInitializer(fixed_image, 
                                                        moving_image, 
                                                        sitk.Euler3DTransform(), 
                                                        sitk.CenteredTransformInitializerFilter.MOMENTS)


    # resample the moving image so that it has the same dimensions and pixel spacing as the fixed image
    resampled_moving_image = sitk.Resample(moving_image, fixed_image, initial_transform)


    # set demons filter and parameters 
    demons_filter = sitk.DiffeomorphicDemonsRegistrationFilter()
    demons_filter.SetNumberOfIterations(iterations) 
    demons_filter.SetSmoothUpdateField(True)
    demons_filter.SetUpdateFieldStandardDeviations(update_field) 
    demons_filter.SetSmoothDisplacementField(True)
    demons_filter.SetStandardDeviations(smooth_displacement)

    # adding callback to track metrics 
    demons_filter.AddCommand(sitk.sitkIterationEvent, lambda: iteration_callback(demons_filter))

    # performing the multiscale registration
    displacement = multiscale_demons(registration_algorithm=demons_filter, 
                        fixed_image = fixed_image, 
                        moving_image = resampled_moving_image,
                        initial_transform = None,
                        shrink_factors = [6,4,2],
                        smoothing_sigmas = [1,1,1])

    # adding the inital transformation and the demons registration together
    displacement_field = sitk.Add(
            displacement,
            sitk.TransformToDisplacementField(
                initial_transform,
                displacement.GetPixelID(),
                displacement.GetSize(),
                displacement.GetOrigin(),
                displacement.GetSpacing(),
                displacement.GetDirection()
            )
        )

    # turning the resulting displacement field into a transformation 
    transformation = sitk.DisplacementFieldTransform(displacement_field)

    # apply transform to the mask
    transformed_mask = sitk.Resample(moving_mask, fixed_image, transformation, sitk.sitkNearestNeighbor, 0.0, moving_mask.GetPixelID())


    # plotting the metric at every iteration to see if it converged (optional)
    if output_graph:
        plt.figure()
        plt.plot(metrics, "k")
        plt.xlabel('Iterations')
        plt.ylabel('Metric')
        plt.grid()
        plt.savefig(output_graph)

    # saving the transformed image
    sitk.WriteImage(transformed_mask, output_image)

def main():
    # Setup description
    description = '''
Utility to perform SimpleITK Demon's deformable registration on two images given a raw atlas image and mask. A new mask will be created for the new image based on the 
given atlas mask. If it does not work sucessfully the first time, the iteration number, update field standard deviations and smooth displacement field standard deviations
can be played with for better results. 
'''
    epilog = '''
Example calls: 
ogoDemonsRegistration raw_atlas.nii.gz atlas_mask.nii.gz new_image.nii.gz output_mask.nii.gz 
ogoDemonsRegistration raw_atlas.nii.gz atlas_mask.nii.gz new_image.nii.gz output_mask.nii.gz --output_graph filename.png --iterations 100
ogoDemonsRegistration raw_atlas.nii.gz atlas_mask.nii.gz new_image.nii.gz output_mask.nii.gz --iterations 100 --update_field 2 --smooth_displacement 2
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoDemonsRegistration",
        description=description,
        epilog=epilog
    )
    parser.add_argument('atlas_path', help = "Full path to where atlas is located, raw atlas image is ARTININ_0009.")
    parser.add_argument('atlas_mask_path', help="Full path to where atlas  mask is located. Pedicles mask is located in ogo/dat.")
    parser.add_argument('new_image', help="Full path to the desired image is, i.e. image without a mask.")
    parser.add_argument('output_image', help="Full path to where transformed mask should go")
    parser.add_argument('-g', '--output_graph', help="Full path to where metrics graph should go")
    parser.add_argument('-i', '--iterations', type=int, default=100, help='Number of iterations to be performed during registration.')
    parser.add_argument('-uf', '--update_field', type=int, default=1, help="Set demons update field standard deviation value, default is 1")
    parser.add_argument('-sd', '--smooth_displacement', type=int, default=1, help="Set demons smooth displacement field standard deviations, default is 1.")


    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('DemonsRegistration', vars(args)))

    # Run program
    DemonsRegistration(**vars(args))


if __name__ == '__main__':
    main()



