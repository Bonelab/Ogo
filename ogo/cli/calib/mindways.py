'''Command line interface for Mindways phantom calibration of CT data.'''

import click
import SimpleITK as sitk
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_csv import write_csv
import numpy as np
from ogo.calibration.mindways_calibration import MindwaysCalibration
import ast


class PythonLiteralOption(click.Option):
    '''A hack for getting variable length options.

    See https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option'''
    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)

@click.command()
@click.argument('ct_file_name', type=click.Path(exists=True))
@click.argument('rods_file_name', type=click.Path(exists=True))
@click.argument('density_file_name', type=click.Path(file_okay=True, writable=True)))
@click.argument('csv_file_name', type=click.Path())
@click.option('--water', nargs=1, cls=PythonLiteralOption, required=True,
    default="[923.2, 1119.52, 1103.57, 1056.95, 1012.25]")
@click.option('--densities', nargs=1, cls=PythonLiteralOption, required=True,
    default="[375.83, 157.05, 58.88, -53.40, -51.83]")
def mindways(ct_file_name, rods_file_name, density_file_name, csv_file_name,
             water, densities):
    '''Determine the calibration equation for a Mindways phantom.
    
    Water and K2HPO4 densities are taken from the QC documentation
    accompanying your specific phantom.'''
    click.echo(echo_arguments('Mindways Calibration', locals()))

    # Read
    click.echo('Reading input {}'.format(ct_file_name))
    ct = sitk.ReadImage(ct_file_name)

    click.echo('Reading input {}'.format(rods_file_name))
    rods = sitk.ReadImage(rods_file_name)

    # Determine number of components
    click.echo('Determining mean intensities')
    filt = sitk.LabelStatisticsImageFilter()
    filt.Execute(ct, rods)

    n_labels = len(filt.GetLabels())
    click.echo('Found {} labels'.format(n_labels))

    if len(water) != len(densities) \
      or len(water) != n_labels-1:
        click.fail('Number of water, K2HPO4, and segmented labels '
        'are not the same ({}, {}, and {})'.format(len(water), len(densities), n_labels-1))

    if len(water) < 2:
        click.fail('Need at least 2 samples to fit a function. Only found {}'.format(len(water)))

    # Get the HU samples
    HU = []
    for label in filt.GetLabels():
        # Skip background label
        if label == 0:
            continue
        
        HU.append(filt.GetMean(label))

    click.echo('Fitting parameters:')
    click.echo('  Water:  {}'.format(water))
    click.echo('  K2HPO4: {}'.format(densities))
    click.echo('  HU:     {}'.format(HU))

    # Perform calibration
    calibrator = MindwaysCalibration()
    calibrator.fit(HU, densities, water)

    click.echo('Found fit:')
    click.echo('  Slope:     {}'.format(calibrator.slope))
    click.echo('  Intercept: {}'.format(calibrator.intercept))
    click.echo('  R^2:       {}'.format(calibrator.r_value**2))

    header = ['CT File Name', 'Rods File Name', 'sigma_ref', 'beta_ref', 'sigma_ct', 'beta_ct', 'slope', 'intercept', 'r_value', 'p_value', 'std_err']
    data = {
        'CT File Name':         ct_file_name,
        'Rods File Name':       rods_file_name,
        'Density File Name':    density_file_name,
        'sigma_ref':            round(calibrator.sigma_ref, 16),
        'beta_ref':             round(calibrator.beta_ref, 16),
        'sigma_ct':             round(calibrator.sigma_ct, 16),
        'beta_ct':              round(calibrator.beta_ct, 16),
        'slope':                round(calibrator.slope, 16),
        'intercept':            round(calibrator.intercept, 16),
        'r_value':              round(calibrator.r_value, 16),
        'p_value':              round(calibrator.p_value, 16),
        'std_err':              round(calibrator.std_err, 16)
    }

    click.echo('Writing complete results to {}'.format(csv_file_name))
    write_csv(data, csv_file_name, header)

    click.echo('Calibrating density file')
    density = calibrator.slope * sitk.Cast(ct, sitkFloat32) + calibrator.intercept

    click.echo('Writing density file to {}'.format(density_file_name))
    sitk.ReadImage(density, density_file_name)