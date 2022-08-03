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
import ogo.cli.Helper as ogo
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
Prints the labels used for ITK-Snap for identifying bones, tissues, and 
calibration rods.

Optionally prints a 128 or 16 palette colormap.
'''
    epilog='''
USAGE: 
ogoPrintLabels --print labels >> ogo-master-labels.txt
ogoPrintLabels --print colors16
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
