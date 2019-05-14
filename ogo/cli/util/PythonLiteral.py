'''Allow Python literals at the command line interface'''

import click
import ast


class PythonLiteral(click.Option):
    '''A hack for getting variable length options.


    This method comes from the following Stack Overflow question:
    https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option'''  # noqa: E501
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:  # noqa: E722
            raise click.BadParameter(value)
