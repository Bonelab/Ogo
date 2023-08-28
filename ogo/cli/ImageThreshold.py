# /------------------------------------------------------------------------------+
# | 19-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import SimpleITK as sitk
import math
import numpy as np
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def ImageThreshold(input_image, output_image, lower, upper, inside, outside, no_output, overwrite):
        
    # Check if output exists and should overwrite
    if os.path.isfile(output_image) and not overwrite and not no_output:
        result = input('File \"{}\" already exists. Overwrite? [y/n]: '.format(output_image))
        if result.lower() not in ['y', 'yes']:
            print('Not overwriting. Exiting...')
            os.sys.exit()
            
    # Read input
    if not os.path.isfile(input_image):
        os.sys.exit('[ERROR] Cannot find file \"{}\"'.format(input_image))

    ogo.message('Reading input CT image to be cropped:')
    ogo.message('      \"{}\"'.format(input_image))
    ct = sitk.ReadImage(input_image)

    # Checking image type
    ogo.message('Input image type is {}'.format(ct.GetPixelIDTypeAsString()))
    image_type = ct.GetPixelIDValue()
    if image_type is sitk.sitkInt8:
        type_min = -128
        type_max = 127
    elif image_type is sitk.sitkUInt8:
        type_min = 0
        type_max = 255
    elif image_type is sitk.sitkInt16:
        type_min = -32768
        type_max = 32767
    elif image_type is sitk.sitkUInt16:
        type_min = 0
        type_max = 65535
    elif image_type is sitk.sitkInt32:
        type_min = -2147483648
        type_max = 2147483647
    elif image_type is sitk.sitkUInt32:
        type_min = 0
        type_max = 4294967295
    else:
        ogo.message('[ERROR] Unexpected image input type.')
    if lower < type_min:
        lower = type_min
        ogo.message('[WARNING] Invalid lower threshold. Changed to {}'.format(lower))
    if upper > type_max:
        upper = type_max
        ogo.message('[WARNING] Invalid upper threshold. Changed to {}'.format(upper))
    
    # Image statistics
    stats = sitk.StatisticsImageFilter()
    stats.Execute(ct)
    
    if lower < stats.GetMinimum():
        ogo.message('Lower ')
    # Report information about input image
    dim = ct.GetSize()
    spacing = ct.GetSpacing()
    origin = ct.GetOrigin()
    phys_dim = [x * y for x, y in zip(dim, spacing)]
    position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    
    guard = '!-------------------------------------------------------------------------------'
    print(guard)
    print('!> dim                            {:>8}  {:>8}  {:>8}'.format(*dim))
    print('!> off                            {:>8}  {:>8}  {:>8}'.format('-', '-', '-'))
    print('!> pos                            {:>8}  {:>8}  {:>8}'.format(*position))
    print('!> element size in mm             {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*spacing))
    print('!> phys dim in mm                 {:>8.4f}  {:>8.4f}  {:>8.4f}'.format(*phys_dim))
    print('!> ')
    print('!> min                            {:>8.0f}'.format(stats.GetMinimum()))
    print('!> max                            {:>8.0f}'.format(stats.GetMaximum()))
    print(guard)
    
    # Threshold
    print(guard)
    print('!> lower                {:>8}'.format(lower))
    print('!> upper                {:>8}'.format(upper))
    print('!> inside                   {:>8}'.format(inside))
    print('!> outside                  {:>8}'.format(outside))
    print(guard)
    
    if lower > upper:
        os.sys.exit('[ERROR] Lower threshold cannot be greater than upper threshold.')
    if inside == outside:
        os.message('[WARNING] Setting inside equal to outside makes no sense!')
        
    ct_thres = sitk.BinaryThresholdImageFilter()
    ct_thres.SetInsideValue(inside)
    ct_thres.SetOutsideValue(outside)
    ct_thres.SetLowerThreshold(lower)
    ct_thres.SetUpperThreshold(upper)
    ct_out = ct_thres.Execute(ct)
    
    # Report stats on resulting image    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    stats.Execute(ct_out,ct)
    
    n_labels = stats.GetNumberOfLabels()
    basename = os.path.basename(input_image)
    
    line_hdr =   'line_{},{},'.format('hdr','name')
    line_units = 'line_{},{},'.format('units','[text]')
    line_data =  'line_{},{},'.format('data',basename)
    line_hdr += '{},{},{},{},{},{},{}'.format('desc','label','total_vol','n_voxels','cx','cy','cz')
    line_units += '{},{},{},{},{},{},{}'.format('[text]','[#]','[mm3]','[N]','[i]','[j]','[k]')
    
    report = ''
    report += '  {:>20s}\n'.format('_______________________________________________________________________Output')    
    report += '  {:>20s} {:>10s} {:>10s} {:>19s}\n'.format('Label','Volume','N Voxels','Centroid')

    if n_labels > 0:
        for label in stats.GetLabels():
            try:
                desc = lb.labels_dict[label]['LABEL']
            except KeyError:
                desc = 'unknown label'
        
            centroid = stats.GetCentroid(label)
            bounding_box = stats.GetBoundingBox(label)
            bb=[0]*3
            bb[0] = bounding_box[0] + int(math.ceil(bounding_box[3]/2))
            bb[1] = bounding_box[1] + int(math.ceil(bounding_box[4]/2))
            bb[2] = bounding_box[2] + int(math.ceil(bounding_box[5]/2))
            line_data += '{},{},{:.1f},{},{},{},{}'.format(desc,label,stats.GetPhysicalSize(label),stats.GetNumberOfPixels(label),bb[0],bb[1],bb[2])
        
            report += '  {:>15s}{:>5s} {:10.1f} {:10d} ({:5d},{:5d},{:5d})\n'\
                      .format('('+desc+')',str(label),stats.GetPhysicalSize(label),stats.GetNumberOfPixels(label),bb[0],bb[1],bb[2])
    else:
        line_data += '{},{},{:.1f},{},{},{},{}'.format('empty','empty',0,0,0,0,0)
        report += '  {:>20s}\n'.format('-- No data --')
        
    print(report)
    print(line_hdr)
    print(line_units)
    print(line_data)
    print()
    
    if not no_output:
        ogo.message('Writing output to file {}'.format(output_image))
        sitk.WriteImage(ct_out, output_image)
    else:
        ogo.message('Writing output to file {} is SUPPRESSED'.format(output_image))
    
    ## Report information about input image
    #dim = ct_out.GetSize()
    #spacing = ct_out.GetSpacing()
    #origin = ct_out.GetOrigin()
    #phys_dim = [x * y for x, y in zip(dim, spacing)]
    #position = [math.floor(x / y) for x, y in zip(origin, spacing)]
    
    ogo.message('Done ogoImageThreshold!')


def main():
    # Setup description
    description = '''
Utility to threshold a NIFTI file.
'''
    epilog = '''
Example calls: 
ogoImageThreshold input.nii.gz --output_image output.nii.gz --upper 1500
ogoImageThreshold input.nii.gz --output_image output.nii.gz --no_output --lower 2500 --upper 5000 --inside 81
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoImageThreshold",
        description=description,
        epilog=epilog
    )
    parser.add_argument('input_image', help='Input image file (*.nii, *.nii.gz)')
    parser.add_argument('output_image', help='Output image file (*.nii, *.nii.gz)')
    parser.add_argument('--lower', type=int, default=250, metavar='VAL', help='Lower threshold (min=-32768, default: %(default)s)')
    parser.add_argument('--upper', type=int, default=3000, metavar='VAL', help='Upper threshold (max=32767, default: %(default)s)')
    parser.add_argument('--inside', type=int, default=127, metavar='VAL', help='Inside value (default: %(default)s)')
    parser.add_argument('--outside', type=int, default=0, metavar='VAL', help='Outside value (default: %(default)s)')
    parser.add_argument('--no_output', action='store_true', help='Suppress writing output (useful if just looking for statistics)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite output without asking')

    # Parse and display
    args = parser.parse_args()
    print(echo_arguments('ImageThreshold', vars(args)))

    # Run program
    ImageThreshold(**vars(args))


if __name__ == '__main__':
    main()
