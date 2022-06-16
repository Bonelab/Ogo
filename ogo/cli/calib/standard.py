'''Command line interface for Mindways phantom calibration of CT data.'''

import click
import SimpleITK as sitk
from collections import OrderedDict
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_txt import write_txt
from ogo.calibration.standard_calibration import StandardCalibration
from ogo.cli.util.PythonLiteral import PythonLiteral
import ogo


@click.command()
@click.argument('ct_file_name', type=click.Path(exists=True))
@click.argument('rods_file_name', type=click.Path(exists=True))
@click.argument('density_file_name',
                type=click.Path(file_okay=True, writable=True))
@click.option('--densities', nargs=1, cls=PythonLiteral, required=True,
              default="[375.83, 157.05, 58.88, -53.40, -51.83]")
@click.option('--calibration_file_name', '-c',
              type=click.Path(file_okay=True, writable=True),
              default='')
def standard(ct_file_name, rods_file_name, density_file_name,
             densities, calibration_file_name):
    '''Determine the calibration equation for a standard phantom.'''
    click.echo(echo_arguments('Standard Calibration', locals()))

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

    if len(densities) < 2:
        click.fail('Need at least 2 samples to fit a function.'
                   ' Only found {}'.format(len(densities)))

    # Get the HU samples
    HU = []
    for label in filt.GetLabels():
        # Skip background label
        if label == 0:
            continue

        HU.append(filt.GetMean(label))

    click.echo('Fitting parameters:')
    click.echo('  Densities: {}'.format(densities))
    click.echo('  HU:        {}'.format(HU))

    # Perform calibration
    calibrator = StandardCalibration()
    calibrator.fit(HU, densities)

    click.echo('Found fit:')
    click.echo('  Slope:     {}'.format(calibrator.slope))
    click.echo('  Intercept: {}'.format(calibrator.intercept))
    click.echo('  R^2:       {}'.format(calibrator.r_value**2))

    click.echo('Calibrating density file')
    density = calibrator.predict(sitk.Cast(ct, sitk.sitkFloat64))
    density = sitk.Cast(density, ct.GetPixelID())

    click.echo('Writing density file to {}'.format(density_file_name))
    sitk.ReadImage(density, density_file_name)

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
