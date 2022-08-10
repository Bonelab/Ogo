'''Command line interface for internal calibration of CT data.'''

import click
import SimpleITK as sitk
import numpy as np
from collections import OrderedDict
from ogo.util.echo_arguments import echo_arguments
from ogo.calibration.internal_calibration import InternalCalibration
from ogo.util.write_txt import write_txt
import ogo


@click.command()
@click.argument('ct_file_name', type=click.Path(exists=True))
@click.argument('internal_file_name', type=click.Path(exists=True))
@click.argument('density_file_name',
                type=click.Path(file_okay=True, writable=True))
@click.option('--calibration_file_name', '-c',
              type=click.Path(file_okay=True, writable=True),
              default='')
def internal(ct_file_name, internal_file_name, density_file_name,
             calibration_file_name):
    '''Perform an internal density calibration

    If you use this work, please cite:
        [1] Michalski et al. 2020 Med Eng Phys.
            https://doi.org/10.1016/j.medengphy.2020.01.009
    '''
    click.echo(echo_arguments('Internal Calibration', locals()))

    labels = {
        'Adipose':          1,
        'Air':              2,
        'Blood':            3,
        'Cortical Bone':    4,
        'Skeletal Muscle':  5
    }

    # Read
    click.echo('Reading input ' + ct_file_name)
    ct = sitk.ReadImage(ct_file_name)

    click.echo('Reading input ' + internal_file_name)
    mask = sitk.ReadImage(internal_file_name)

    # Determine number of components
    click.echo('Determining mean intensities')
    filt = sitk.LabelStatisticsImageFilter()
    filt.Execute(ct, mask)

    n_labels = len(filt.GetLabels())
    click.echo('  Found {} labels'.format(n_labels))

    for label, value in labels.items():
        if not filt.HasLabel(value):
            click.echo(
                '  Could not find values for label \"{}\" ({})'.format(
                    label, value
                ), err=True)
            raise click.Abort()
    click.echo('  Expected labels found')

    click.echo('Computing calibration parameters')
    calib = InternalCalibration(
        adipose_hu=filt.GetMean(labels['Adipose']),
        air_hu=filt.GetMean(labels['Air']),
        blood_hu=filt.GetMean(labels['Blood']),
        bone_hu=filt.GetMean(labels['Cortical Bone']),
        muscle_hu=filt.GetMean(labels['Skeletal Muscle'])
    )
    calib.fit()

    voxel_volume = np.prod(ct.GetSpacing())

    click.echo('Calibrating input file')
    den = calib.predict(sitk.Cast(ct, sitk.sitkFloat64), voxel_volume)
    den = sitk.Cast(den, ct.GetPixelID())

    click.echo('Writing result to ' + density_file_name)
    sitk.WriteImage(den, density_file_name)

    if calibration_file_name:
        click.echo('Saving calibration parameters to file ' +
                   calibration_file_name)
        entry = OrderedDict([
            ('ct_file_name',            ct_file_name),
            ('internal_file_name',      internal_file_name),
            ('density_file_name',       density_file_name),
            ('calibration_file_name',   calibration_file_name),
            ('Ogo Version',             ogo.__version__)
        ])
        entry.update(calib.get_dict())
        write_txt(entry, calibration_file_name)
