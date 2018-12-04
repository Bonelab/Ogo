'''Command line interface for isotropic resampling of QCT data'''

import click
import SimpleITK as sitk
from ogo.util.echo_arguments import echo_arguments
import numpy as np


@click.command()
@click.argument('input_file_name', type=click.Path(exists=True))
@click.argument('output_file_name',
                type=click.Path(file_okay=True, writable=True))
@click.option('-r', '--isotropic_resolution', type=float, default=1.0)
def resample(input_file_name, output_file_name,
             isotropic_resolution):
    '''Resample a CT image to an isotropic resolution'''
    click.echo(echo_arguments('Isotropic Resampling', locals()))

    # Read
    click.echo('Reading input {}'.format(input_file_name))
    image = sitk.ReadImage(input_file_name)

    # Determine sigma
    click.echo('Resampling to a resolution of {}'.format(isotropic_resolution))
    sigma = isotropic_resolution * np.sqrt(2 * np.log(2)) / np.pi
    click.echo('  For a half-power cut-off frequency of {:0.3f}, '
               'a standard deviation of {} is being used.'.format(
                  1/(2.0*isotropic_resolution), sigma))

    # Filter each sampling directions as needed
    for i, spacing in enumerate(image.GetSpacing()):
        if spacing > isotropic_resolution:
            click.echo('  No antialiasing in direction {}'.format(i))
            continue

        click.echo('  Antialiasing in direction {}'.format(i))
        image = sitk.RecursiveGaussian(
          image,
          sigma, False,
          sitk.RecursiveGaussianImageFilter.ZeroOrder,
          i
        )

    # Determine the output size
    resolution = [isotropic_resolution for i in range(image.GetDimension())]
    size = [int(np.ceil(s * i / o)) for o, i, s in
            zip(resolution, image.GetSpacing(), image.GetSize())]
    click.echo('  Input Size:  {}'.format(image.GetSize()))
    click.echo('  Output Size: {}'.format(size))

    click.echo('Finding minimum in image')
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    vox_min = stats.GetMinimum()
    click.echo('  Minimum intensity: {}'.format(vox_min))

    click.echo('Resampling...')
    transform = sitk.Euler3DTransform()
    transform.SetIdentity()

    click.echo('  Using BSpline interpolation')

    output = sitk.Resample(
      image,
      size,
      transform,
      sitk.sitkBSpline,
      image.GetOrigin(),
      resolution,
      image.GetDirection(),
      vox_min,
      image.GetPixelID()
    )

    click.echo('Writing output to file {}'.format(output_file_name))
    sitk.WriteImage(output, output_file_name)
