"""Image I/O helpers for FE sidecar files."""

from pathlib import Path


def write_vtk_image_with_sitk_geometry(vtk_image, output_file):
    """Write a VTK image with spacing/origin preserved for SimpleITK readers."""
    import SimpleITK as sitk
    from vtk.util.numpy_support import vtk_to_numpy

    output_file = Path(output_file)
    dims = vtk_image.GetDimensions()
    data = vtk_to_numpy(vtk_image.GetPointData().GetScalars()).reshape(dims, order="F")
    sitk_image = sitk.GetImageFromArray(data.transpose(2, 1, 0))
    sitk_image.SetSpacing(tuple(float(v) for v in vtk_image.GetSpacing()))
    sitk_image.SetOrigin(tuple(float(v) for v in vtk_image.GetOrigin()))
    sitk_image.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    output_file.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(sitk_image, str(output_file))
    return output_file
