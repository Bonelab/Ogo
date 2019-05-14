'''Command line interface for Mindways phantom calibration of CT data.'''

import click
import SimpleITK as sitk
from ogo.util.echo_arguments import echo_arguments
from ogo.util.write_csv import write_csv
from ogo.calibration.standard_calibration import StandardCalibration
from ogo.cli.util.PythonLiteral import PythonLiteral


@click.command()
@click.argument('ct_file_name', type=click.Path(exists=True))
@click.argument('rods_file_name', type=click.Path(exists=True))
@click.argument('density_file_name',
                type=click.Path(file_okay=True, writable=True))
@click.option('--densities', nargs=1, cls=PythonLiteral, required=True,
              default="[375.83, 157.05, 58.88, -53.40, -51.83]")
@click.argument('csv_file_name', type=click.Path())
def standard(ct_file_name, rods_file_name, density_file_name,
             densities, csv_file_name):
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

    header = [
        'CT File Name', 'Rods File Name', 'Density File Name',
        'slope', 'intercept', 'r_value', 'p_value', 'std_err'
    ]
    data = {
        'CT File Name':         ct_file_name,
        'Rods File Name':       rods_file_name,
        'Density File Name':    density_file_name,
        'slope':                round(calibrator.slope, 16),
        'intercept':            round(calibrator.intercept, 16),
        'r_value':              round(calibrator.r_value, 16),
        'p_value':              round(calibrator.p_value, 16),
        'std_err':              round(calibrator.std_err, 16)
    }

    click.echo('Writing complete results to {}'.format(csv_file_name))
    write_csv(data, csv_file_name, header)

    click.echo('Calibrating density file')
    density = calibrator.slope * sitk.Cast(ct, sitk.sitkFloat32) + \
        calibrator.intercept

    click.echo('Writing density file to {}'.format(density_file_name))
    sitk.ReadImage(density, density_file_name)
