# /------------------------------------------------------------------------------+
# | 06-JUL-2022                                                                  |
# | Copyright (c) Bone Imaging Laboratory                                        |
# | All rights reserved                                                          |
# | bonelab@ucalgary.ca                                                          |
# +------------------------------------------------------------------------------+

'''Utility function for writing to a csv file'''

import os


def write_csv(entry, csv_file, header, delimiter=','):
    '''Write data to a CSV file.

    This function handles creating the CSV file if it doesn''t exist,
    appending to said file, and managing header formatting.

    Note that if you want specific formatting on numerics - i.e. precision
    - the numerical value should be passed as a string and precision set
    before calling this function.'''
    # Create formating strings
    format_header = ['{{{}}}'.format(x) for x in header]
    format_string = delimiter.join(format_header)

    # If we doesn't exist, need to write header
    if not os.path.exists(csv_file):
        with open(csv_file, 'w') as f:
            f.write(delimiter.join(header) + os.linesep)

    # Write entry
    with open(csv_file, 'a') as f:
        f.write(format_string.format(**entry) + os.linesep)
