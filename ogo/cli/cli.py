'''Command line interface entry point'''

import click
from .img.img import img


@click.group()
def main():
    '''Osteoporosis is eluse - are you going to help catch it?'''
    pass


main.add_command(img)
