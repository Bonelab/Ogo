'''Command line interface for calibrating CT data'''

import click
from .mindways import mindways
from .standard import standard
from .internal import internal


@click.group()
def calib():
    '''Command line interface for calibrating CT data'''
    pass


calib.add_command(mindways)
calib.add_command(standard)
calib.add_command(internal)
