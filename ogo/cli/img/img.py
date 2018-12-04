'''Command line interface for basic image processing'''

import click
from .isotropic_resampling import resample


@click.group()
def img():
    '''Various commands for basic image processing'''
    click.echo


img.add_command(resample)
