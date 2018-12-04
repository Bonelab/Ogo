'''Utility function for echoing command line arguments'''

import os


def echo_arguments(title, args):
    '''Echo the arguments passed to a function'''
    # Title
    message = title + ':' + os.linesep

    # Determine padding size
    max_length = 0
    for key in args:
        max_length = max(max_length, len(key))
    formatter = '  {{arg:<{spacing}}}{{value}}'.format(spacing=max_length+2)

    # Print arguments and values
    for key in args:
        message += formatter.format(arg=key, value=args[key]) + os.linesep

    return message
