'''Command line interface for calibrating CT data'''

import click
from .mindways import mindways
from .standard import standard


@click.group()
def calib():
    '''Command line interface for calibrating CT data'''
    pass


calib.add_command(mindways)
calib.add_command(standard)
