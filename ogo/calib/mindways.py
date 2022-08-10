'''Command line interface for Mindways phantom calibration of CT data.'''

import click
import SimpleITK as sitk
from collections import OrderedDict
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_txt import write_txt
from ogo.calibration.mindways_calibration import MindwaysCalibration
from ogo.cli.util.PythonLiteral import PythonLiteral
import ogo


@click.command()
@click.argument('ct_file_name', type=click.Path(exists=True))
@click.argument('rods_file_name', type=click.Path(exists=True))
@click.argument('density_file_name',
                type=click.Path(file_okay=True, writable=True))
@click.option('--calibration_file_name', '-c',
              type=click.Path(file_okay=True, writable=True),
              default='')
@click.option('--water', nargs=1, cls=PythonLiteral, required=True,
              default="[923.2, 1119.52, 1103.57, 1056.95, 1012.25]")
@click.option('--densities', nargs=1, cls=PythonLiteral, required=True,
              default="[375.83, 157.05, 58.88, -53.40, -51.83]")
def mindways(ct_file_name, rods_file_name, density_file_name,
             calibration_file_name,
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

    if len(water) != len(densities) or len(water) != n_labels-1:
        click.fail(
            'Number of water, K2HPO4, and segmented labels are not the same '
            '({}, {}, and {})'.format(len(water), len(densities), n_labels-1))

    if len(water) < 2:
        click.fail('Need at least 2 samples to fit a function.'
                   ' Only found {}'.format(len(water)))

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

    click.echo('Calibrating density file')
    density = calibrator.predict(sitk.Cast(ct, sitk.sitkFloat64))
    density = sitk.Cast(density, ct.GetPixelID())

    click.echo('Writing density file to ' + density_file_name)
    sitk.WriteImage(density, density_file_name)

    if calibration_file_name:
        click.echo('Saving calibration parameters to file ' +
                   calibration_file_name)
        entry = OrderedDict([
            ('ct_file_name',            ct_file_name),
            ('rods_file_name',          rods_file_name),
            ('density_file_name',       density_file_name),
            ('calibration_file_name',   calibration_file_name),
            ('Ogo Version',             ogo.__version__)
        ])
        entry.update(calibrator.get_dict())
        write_txt(entry, calibration_file_name)
