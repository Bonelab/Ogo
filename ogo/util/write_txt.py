# /------------------------------------------------------------------------------+
# | 22-AUG-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

'''Utility function for writing to a text file'''

import os


def write_txt(entry, filename, mode='w', delimiter=','):
    '''Write data to a text file.

    `entry` should be a collections.OrderedDict object to preserve the
    ordering of the variables.

    `mode` can be 'a' for append or 'w' for overwrite.

    Note that if you want specific formatting on numerics - i.e. precision
    - the numerical value should be passed as a string and precision set
    before calling this function.'''
    # Write entry
    with open(filename, mode) as f:
        for key, value in entry.items():
            f.write(delimiter.join([str(key), str(value)]) + os.linesep)
