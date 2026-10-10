"""Small VTK image/NumPy conversion helpers used by FE workflows."""


def vtk_image_to_numpy(vtk_image, *, processing_order=False, copy=False):
    """Return image scalars as a NumPy array.

    By default the returned array is in VTK image dimension order.  Set
    ``processing_order`` when using FE boundary masks that are processed as
    z/y/x NumPy arrays.
    """
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy

    scalars = vtk_image.GetPointData().GetScalars()
    try:
        array = vtk_to_numpy(scalars)
    except ValueError:
        array = np.asarray([scalars.GetTuple1(i) for i in range(scalars.GetNumberOfTuples())])
    array = array.reshape(vtk_image.GetDimensions(), order="F")
    if processing_order:
        array = np.swapaxes(array, 0, 2)
    if copy:
        array = array.copy()
    return array


def numpy_to_vtk_image(array, template_vtk_image, vtk_array_type=None, *, processing_order=False):
    """Create a VTK image from an array and preserve the template image geometry."""
    import numpy as np
    from vtk.util.numpy_support import numpy_to_vtk

    final_array = np.swapaxes(array, 0, 2) if processing_order else array
    if vtk_array_type is None:
        vtk_array_type = template_vtk_image.GetScalarType()
    vtk_array = numpy_to_vtk(final_array.ravel(order="F"), deep=True, array_type=vtk_array_type)

    output = template_vtk_image.NewInstance()
    output.DeepCopy(template_vtk_image)
    output.SetExtent(template_vtk_image.GetExtent())
    output.SetOrigin(template_vtk_image.GetOrigin())
    output.SetSpacing(template_vtk_image.GetSpacing())
    output.GetPointData().SetScalars(vtk_array)
    return output
