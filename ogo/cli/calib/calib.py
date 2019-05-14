'''Command line interface for calibrating CT data'''

import click
from .mindways import mindways


@click.group()
def calib():
    '''Command line interface for calibrating CT data'''
    pass


calib.add_command(mindways)
