# /------------------------------------------------------------------------------+
# | 2-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

# Imports
import argparse
import os
import vtk
from ogo.util.echo_arguments import echo_arguments
import ogo.util.Helper as ogo
import ogo.dat.OgoMasterLabels as lb

def PrintLabels(echo):

    if echo in 'labels':
        print(lb.labels_hdr)
        print("\n".join("{:5d} {:4d} {:4d} {:4d}  {:.1f} {:4d} {:4d}  \"{}\"".format(k,v['RGB'][0],\
                                                                                       v['RGB'][1],\
                                                                                       v['RGB'][2],\
                                                                                       v['A'],\
                                                                                       v['VIS'],\
                                                                                       v['IDX'],\
                                                                                       v['LABEL'])\
                                                           for k, v in lb.labels_dict.items()))
    
    elif echo in 'colors128':
        print(lb.colors128_hdr)
        print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],\
                                                                    v['RGB'][1],\
                                                                    v['RGB'][2],\
                                                                    v['DESC'])\
                                                           for k, v in lb.colors128_dict.items()))
        
    elif echo in 'colors16':
        print(lb.colors16_hdr)
        print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],\
                                                                    v['RGB'][1],\
                                                                    v['RGB'][2],\
                                                                    v['DESC'])\
                                                           for k, v in lb.colors16_dict.items()))
        
    else:
        os.sys.exit('[ERROR] Unknown argument: \"{}\"'.format(echo))
        
    # labels
    #print(labels_hdr)
    #print("\n".join("{:5d} {:4d} {:4d} {:4d}  {:.1f} {:4d} {:4d}  \"{}\"".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['A'],v['VIS'],v['IDX'],v['LABEL']) \
    #                                                   for k, v in labels_dict.items()))
    
    # color16
    #print(colors16_hdr)
    #print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['DESC']) \
    #                                                   for k, v in colors16_dict.items()))
    
    # color128
    #print(colors128_hdr)
    #print("\n".join("{:5d} {:4d}  {:4d}  {:4d}   {:s}".format(k,v['RGB'][0],v['RGB'][1],v['RGB'][2],v['DESC']) \
    #                                                   for k, v in colors128_dict.items()))
    
        
def main():
    # Setup description
    description='''
Prints labels that are the standard for all Ogo applications to identify 
bones, tissues, and calibration rodes. 

The master list of labels is kept in OgoMasterLabels.py. Any changes or 
additions should be made in that document. This application reads that 
document to produce output. 

A text file created by the output from this application can be read in by 
ITK-Snap:

[Segmentation] --> [Import Label Descriptions]

This application also contains information about 16 and 128 colour maps 
that can be useful if/when updating OgoMasterLabels.py.

'''
    epilog='''
Example calls: 
ogoPrintLabels --echo labels >> ogo-master-labels.txt
ogoPrintLabels --echo colors16
'''

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoPrintLabels",
        description=description,
        epilog=epilog
    )
    parser.add_argument('--echo', default='labels', choices=['labels','colors16','colors128'], help='Specify output (default: %(default)s)')

    print()

    # Parse and display
    args = parser.parse_args()
    #print(echo_arguments('PrintLabels', vars(args)))

    # Run program
    PrintLabels(**vars(args))

if __name__ == '__main__':
    main()
