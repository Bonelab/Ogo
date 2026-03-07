#####
# ogo_sideways_fall_fe.py
#
# This script sets up the sideways fall FE model on the hip from the density (K2HPO4)
# calibrated image. This script sets up the model for either a left or right femur, as
# specified by the user. The analysis resamples the image to isotropic voxels, transforms
# the image, applies the bone mask and bins the data. It then creates the FE model for
# solving using FAIM (>v8.0, Numerics Solutions Ltd, Calgary, Canada - Steven Boyd).
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 29, 2019
# Modified to Py3: March 25, 2020
# Patched for material-table API compatibility and pre-transform padding
#####

script_version = 1.1

##
# Import the required modules
import ogo.util.Helper as ogo
import os
import sys
import argparse
import vtk
import vtkbone
import numpy as np

from ogo.util.echo_arguments import echo_arguments


def remove_extension(filename):
    while True:
        filename, ext = os.path.splitext(filename)
        if not ext:
            break
    return filename


def pad_image(image, pad_vox=20, pad_value=0):
    """
    Pad vtkImageData by a fixed number of voxels on all sides.

    Parameters
    ----------
    image : vtk.vtkImageData
        Input image.
    pad_vox : int
        Number of voxels to pad on each side in x, y, z.
    pad_value : numeric
        Constant value used for padding.

    Returns
    -------
    vtk.vtkImageData
        Padded image with updated origin and same spacing.
    """
    extent = list(image.GetExtent())  # [xmin, xmax, ymin, ymax, zmin, zmax]
    spacing = image.GetSpacing()
    origin = list(image.GetOrigin())

    new_extent = [
        extent[0] - pad_vox, extent[1] + pad_vox,
        extent[2] - pad_vox, extent[3] + pad_vox,
        extent[4] - pad_vox, extent[5] + pad_vox,
    ]

    new_origin = [
        origin[0] - pad_vox * spacing[0],
        origin[1] - pad_vox * spacing[1],
        origin[2] - pad_vox * spacing[2],
    ]

    padded = vtk.vtkImageConstantPad()
    padded.SetInputData(image)
    padded.SetOutputWholeExtent(*new_extent)
    padded.SetConstant(pad_value)
    padded.Update()

    out = vtk.vtkImageData()
    out.DeepCopy(padded.GetOutput())
    out.SetOrigin(new_origin)
    out.SetSpacing(spacing)

    # Preserve direction matrix if available
    if hasattr(image, "GetDirectionMatrix") and hasattr(out, "SetDirectionMatrix"):
        try:
            out.SetDirectionMatrix(image.GetDirectionMatrix())
        except Exception:
            pass

    return out


##
# Start script
def sidewaysFallFe(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask

    mask_threshold = args.mask_threshold
    iso_resolution = args.iso_resolution
    femur_side = args.femur_side
    elastic_exponent = args.elastic_exponent
    elastic_Emax = args.elastic_Emax
    poissons_ratio = args.poissons_ratio
    pmma_E = args.pmma_E
    pmma_v = args.pmma_v
    pmma_thick = args.pmma_thick
    pmma_mat_id = args.pmma_mat_id
    fe_displacement = args.fe_displacement
    left_femur_reference = args.left_femur_reference
    right_femur_reference = args.right_femur_reference
    output_file = args.output_file
    bone_yield_compression = args.bone_yield_compression
    bone_yield_tension = args.bone_yield_tension
    pmma_yield_compression = args.pmma_yield_compression
    pmma_yield_tension = args.pmma_yield_tension
    pad_vox = args.pad_vox

    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)
    script_name = sys.argv[0]

    if femur_side == 1:
        N88_fileName = output_file.replace('.n88model', "_LT_FEMUR_SF.n88model")
        rot_z = 90
    elif femur_side == 2:
        N88_fileName = output_file.replace('.n88model', "_RT_FEMUR_SF.n88model")
        rot_z = -90
    else:
        ogo.message("Femur side not recognized. Terminating...")
        sys.exit()

    ##
    # Message the input parameters to the terminal
    ogo.message("Image Path: %s" % image_pathname)
    ogo.message("Image File: %s" % image_basename)
    ogo.message("Mask Path: %s" % mask_pathname)
    ogo.message("Mask File: %s" % mask_basename)
    ogo.message("Isotropic Voxel Size: %8.4f" % iso_resolution)
    ogo.message("Pre-transform Padding [vox]: %d" % pad_vox)
    if femur_side == 1:
        ogo.message("Femur Side for Model: Left")
    elif femur_side == 2:
        ogo.message("Femur Side for Model: Right")
    ogo.message("Bone Power Law Exponent: %8.4f" % elastic_exponent)
    ogo.message("Bone Power Law Constant [MPa]: %8.4f" % elastic_Emax)
    ogo.message("Bone Poissons Ratio: %1.1f" % poissons_ratio)
    ogo.message("PMMA Elastic Modulus [MPa]: %8.4f" % pmma_E)
    ogo.message("PMMA Poissons Ration: %1.1f" % pmma_v)
    ogo.message("PMMA Thickness [mm]: %8.4f" % pmma_thick)
    ogo.message("PMMA Material ID: %d" % pmma_mat_id)
    ogo.message("Applied Displacement [mm]: %8.4f" % fe_displacement)

    ##
    # Read input image
    ogo.message("Reading calibrated image...")
    imageData = ogo.readNii(image)

    ##
    # Read bone mask
    ogo.message("Reading bone mask...")
    maskData = ogo.readNii(mask)
    maskThres = ogo.maskThreshold(maskData, mask_threshold)

    ##
    # Resampling the images to isotropic voxel size
    ogo.message("Resampling the input image and bone mask to isotropic...")
    image_resample = ogo.imageResample(imageData, iso_resolution)
    mask_resample = ogo.imageResample(maskThres, iso_resolution)

    ##
    # Pad image and mask before rotation / ICP to reduce clipping
    ogo.message("Padding image and bone mask before transformation...")
    image_pad = pad_image(image_resample, pad_vox=pad_vox, pad_value=0)
    mask_pad = pad_image(mask_resample, pad_vox=pad_vox, pad_value=0)

    ##
    # Pre-rotate the image and mask for better alignment
    ogo.message("Pre-rotating image and mask...")
    image_rot, mask_rot = ogo.preRotateImage(image_pad, mask_pad, rot_z)

    ##
    # Align the input femur with the reference model
    ogo.message("Aligning input with reference model...")
    mask_surface = ogo.marchingCubes(mask_rot)
    if femur_side == 1:
        ref_poly = ogo.readPolyData(left_femur_reference)
    elif femur_side == 2:
        ref_poly = ogo.readPolyData(right_femur_reference)
    else:
        print("Error: Femur Side not defined. Terminating...")
        sys.exit()

    icp = ogo.iterativeClosestPoint(ref_poly, mask_surface)

    ogo.message("Applying the transformation to the image and mask...")
    image_trans = ogo.applyTransform(image_rot, icp)
    mask_trans = ogo.applyTransform(mask_rot, icp)

    # replace any negative values less than -31 to be equivalent to -31.
    # -31 is used as that converts to a minimum elastic modulus value of 0.1 MPa
    # K2HPO4 den = -31 mg/cc => Ash den = 6 mg/cc => E = 0.1 MPa
    image_thres = ogo.bmd_preprocess(image_trans, -31)
    image_ash = ogo.bmd_K2hpo4ToAsh(image_thres)

    binned_image, bin_centers = ogo.density2materialID(image_ash, n_bins=128)

    ##
    # Set up the FE model.
    ogo.message("Setting up the Finite Element Model...")
    # Cast the image to Short to round float values to nearest whole number
    cast_image = ogo.cast2short(binned_image)
    cast_mask = ogo.cast2unsignchar(mask_trans)

    # Apply the mask to the bone
    ogo.message("Applying the bone mask to the image...")
    bone_image = ogo.applyMask(cast_image, cast_mask)

    # Ensure connectivity
    ogo.message("Performing image connectivity...")
    conn = ogo.imageConnectivity(bone_image)

    change = ogo.changeInfo(conn)

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh = ogo.Image2Mesh(change)

    ##
    # Set up the Material Table
    ogo.message("Setting up the Finite Element Material Table...")
    material_table = vtkbone.vtkboneMaterialTable()
    material_table = ogo.add_bone_material_depreciated(
        material_table,
        bin_centers,
        elastic_Emax=elastic_Emax,
        elastic_exponent=elastic_exponent,
        mu=poissons_ratio,
        bone_yield_compression=bone_yield_compression,
        bone_yield_tension=bone_yield_tension,
        bin_range=(0, len(bin_centers)),
        material_name='BoneMat'
    )
    material_table = ogo.add_pmma_material(
        material_table,
        pmma_mat_id,
        pmma_E,
        pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression
    )

    ##
    # Create the preliminary FE model
    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()
    model_bounds = model.GetBounds()
    ogo.message("Model Bounds: %s" % str(model_bounds))

    # Define support vectors for later application
    top_support_vector = (0, 1, 0)
    bottom_support_vector = (0, -1, 0)
    side_support_vector = (0, 0, -1)

    # Define the rotation angle in degrees
    theta = np.radians(-20)

    # Define the rotation matrix around the x-axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta),  np.cos(theta)]
    ])

    # Apply the rotation
    side_support_vector = Rx @ side_support_vector

    ogo.message("Determining Femoral Head ROI...")
    femoral_head_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[2]),
        float(model_bounds[2] + (model_bounds[3] - model_bounds[2]) / 5),
        float(model_bounds[4] + (model_bounds[5] - model_bounds[4]) / 2),
        float(model_bounds[5])
    )
    ogo.message("Extracted Femoral Head Bounds: %s" % str(femoral_head_bounds))
    femoral_head_model = ogo.extractBox(femoral_head_bounds, model)
    femoral_head_model_bounds = femoral_head_model.GetBounds()
    ogo.message("Extracted Femoral Head Model Bounds: %s" % str(femoral_head_model_bounds))
    femoral_head_model_bounds = (
        femoral_head_model_bounds[0],
        femoral_head_model_bounds[1] - 1,
        femoral_head_model_bounds[2],
        femoral_head_model_bounds[3] - 1,
        femoral_head_model_bounds[4],
        femoral_head_model_bounds[5] - 1
    )

    ##
    # Extract Greater Trochanter ROI
    ogo.message("Determining Greater Trochanter ROI...")
    greater_trochanter_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[3] - (model_bounds[3] - model_bounds[2]) / 20),
        float(model_bounds[3]),
        float(model_bounds[4]),
        float(model_bounds[5])
    )
    ogo.message("Extracted Greater Trochanter Bounds: %s" % str(greater_trochanter_bounds))
    greater_trochanter_model = ogo.extractBox(greater_trochanter_bounds, model)
    greater_trochanter_model_bounds = greater_trochanter_model.GetBounds()
    ogo.message("Extracted Greater Trochanter Model Bounds: %s" % str(greater_trochanter_model_bounds))

    greater_trochanter_model_bounds = (
        greater_trochanter_model_bounds[0],
        greater_trochanter_model_bounds[1] - 1,
        greater_trochanter_model_bounds[2],
        greater_trochanter_model_bounds[3] - 1,
        greater_trochanter_model_bounds[4],
        greater_trochanter_model_bounds[5] - 1
    )

    ##
    # Define values for adding PMMA blocks
    inval = pmma_mat_id
    outval = 0

    # Thickness to be added to PMMA cap.
    origin = image_resample.GetOrigin()
    spacing = image_resample.GetSpacing()
    thickness = int(pmma_thick / spacing[1])

    ##
    # Creates the femoral head PMMA cap
    ogo.message("Creating Femoral Head PMMA Cap...")
    femoralHeadPMMA = ogo.femoralHeadPMMA(
        femoral_head_model_bounds,
        spacing,
        origin,
        inval,
        outval,
        thickness,
        pmma_mat_id
    )

    ##
    # Creates the greater trochanter PMMA cap
    ogo.message("Creating Greater Trochanter PMMA Cap...")
    greaterTrochanterPMMA = ogo.greaterTrochanterPMMA(
        greater_trochanter_model_bounds,
        spacing,
        origin,
        inval,
        outval,
        thickness,
        pmma_mat_id
    )

    ##
    # Combines the PMMA cap images with the model image
    ogo.message("Combine PMMA Cap Images with Model Image...")
    combinedImage = ogo.combineImageData_SF(conn, femoralHeadPMMA, greaterTrochanterPMMA, pmma_mat_id)

    ##
    # Mesh the final image and create Finite Element Model
    ogo.message("Performing Image connectivity...")
    conn2 = ogo.imageConnectivity(combinedImage)

    ogo.message("Meshing image data to elements...")
    mesh2 = ogo.Image2Mesh(conn2)

    ##
    # Create the final FE model
    ogo.message("Constructing the Finite Element Model...")
    model2 = ogo.applyTestBase(mesh2, material_table)
    model2.ComputeBounds()
    model2_bounds = model2.GetBounds()
    ogo.message("Model 2 Bounds: %s" % str(model2_bounds))

    ##
    # Determine Femoral Head PMMA Cap support nodes
    ogo.message("Determining Femoral Head PMMA Cap nodes...")
    fh_pmma_bounds = (
        model2_bounds[0],
        model2_bounds[1],
        model2_bounds[2],
        model2_bounds[2] + 1,
        model2_bounds[4],
        model2_bounds[5]
    )
    ogo.message("Femoral Head PMMA Cap Bounds: %s" % str(fh_pmma_bounds))
    fh_pmma_model = ogo.extractBox(fh_pmma_bounds, model2)
    fh_pmma_model_bounds = fh_pmma_model.GetBounds()
    ogo.message("Femoral Head PMMA Cap Model Bounds: %s" % str(fh_pmma_model_bounds))

    fh_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        fh_pmma_visible_node_IDS,
        fh_pmma_model,
        bottom_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Femoral Head PMMA Cap."
                % fh_pmma_visible_node_IDS.GetNumberOfTuples())
    fh_pmma_visible_node_IDS.SetName("Femoral_Head_PMMA_Nodes")
    model2.AddNodeSet(fh_pmma_visible_node_IDS)

    ##
    # Determine Greater Trochanter PMMA Cap support nodes
    ogo.message("Determining Greater Trochanter PMMA Cap support nodes...")
    gt_pmma_bounds = (
        model2_bounds[0],
        model2_bounds[1],
        model2_bounds[3] - 1,
        model2_bounds[3],
        model2_bounds[4],
        model2_bounds[5]
    )
    ogo.message("Greater Trochanter PMMA Cap Bounds: %s" % str(gt_pmma_bounds))
    gt_pmma_model = ogo.extractBox(gt_pmma_bounds, model2)
    gt_pmma_model_bounds = gt_pmma_model.GetBounds()
    ogo.message("Greater Trochanter PMMA Cap Model Bounds: %s" % str(gt_pmma_model_bounds))

    gt_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        gt_pmma_visible_node_IDS,
        gt_pmma_model,
        top_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Greater Trochanter PMMA Cap."
                % gt_pmma_visible_node_IDS.GetNumberOfTuples())

    gt_pmma_visible_node_IDS.SetName("Greater_Trochanter_PMMA_Nodes")
    model2.AddNodeSet(gt_pmma_visible_node_IDS)

    ##
    # Distal femur support nodes
    ogo.message("Determining distal femur nodes...")
    distal_femur_bounds = (
        model2_bounds[0],
        model2_bounds[1],
        model2_bounds[2],
        model2_bounds[3],
        model2_bounds[4],
        model2_bounds[4] + 10
    )
    ogo.message("Distal Femur Bounds: %s" % str(distal_femur_bounds))
    df_model = ogo.extractBox(distal_femur_bounds, model2)
    df_model_bounds = df_model.GetBounds()
    ogo.message("Distal Femur Model Bounds: %s" % str(df_model_bounds))

    df_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        df_visible_node_IDS,
        df_model,
        side_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Distal Femur."
                % df_visible_node_IDS.GetNumberOfTuples())

    df_visible_node_IDS.SetName("Distal_Femur_Nodes")
    model2.AddNodeSet(df_visible_node_IDS)

    ##
    # Apply Boundary conditions to PMMA caps at specific sites
    ogo.message("Applying boundary conditions to Femoral Head PMMA cap...")
    model2.ApplyBoundaryCondition(
        "Femoral_Head_PMMA_Nodes",
        vtkbone.vtkboneConstraint.SENSE_Y,
        -fe_displacement,
        "top_displacement"
    )

    dirs = [
        vtkbone.vtkboneConstraint.SENSE_Y,
        vtkbone.vtkboneConstraint.SENSE_X,
        vtkbone.vtkboneConstraint.SENSE_Z
    ]
    labels = ['y', 'x', 'z']

    for SENSE, label in zip(dirs, labels):
        ogo.message("Applying boundary conditions to Greater Trochanter PMMA cap...")
        model2.ApplyBoundaryCondition(
            "Greater_Trochanter_PMMA_Nodes",
            SENSE,
            0,
            f"bottom_fixed_{label}_PMMA"
        )

        ogo.message("Applying boundary conditions to Distal Femur...")
        model2.ApplyBoundaryCondition(
            "Distal_Femur_Nodes",
            SENSE,
            0,
            f"bottom_fixed_{label}"
        )

    ##
    # Post Processing parameters
    ogo.message("Setting up Post Processing Parameters...")
    info = model2.GetInformation()
    pp_node_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_NODE_SETS()
    pp_elem_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_ELEMENT_SETS()

    for setname in ["Femoral_Head_PMMA_Nodes", "Greater_Trochanter_PMMA_Nodes", "Distal_Femur_Nodes"]:
        pp_node_sets_key.Append(info, setname)
        elementSet = model2.GetAssociatedElementsFromNodeSet(setname)
        model2.AddElementSet(elementSet)
        pp_elem_sets_key.Append(info, setname)

    model2.AppendHistory("Created by %s version %s." % (script_name, script_version))

    ##
    # Write out n88model file
    ogo.message("Writing out n88model file: %s" % N88_fileName)
    ogo.writeN88Model(model2, N88_fileName, image_pathname)

    ogo.message("Done writing n88model.")
    ogo.message("End of Script.")
    sys.exit()


def main():
    description = '''
This script sets up the sideways fall FE model on the hip from the
density (K2HPO4) calibrated image. This script sets up the model for either a left or
right femur, as specified by the user. The analysis resamples the image to isotropic
voxels, transforms the image, applies the bone mask and bins the data. It then
creates the FE model for solving using FAIM (v8.1, Numerics Solutions Ltd, Calgary,
Canada - Steven Boyd).

Input: Calibrated K2HPO4 Image (*.nii), Bone Mask (*_MASK.nii)

Optional Parameters:
1) Isotropic resample voxel size
2) Left or Right Femur Bone Mask

Output: N88 Model (*.n88model)
'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="OgoSidewaysFallFe",
        description=description
    )

    parser.add_argument("calibrated_image",
                        help="*_K2HPO4.nii image file")
    parser.add_argument("bone_mask",
                        help="*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type=int, default=1,
                        help="Set the threshold value to extract the bone of interest from the mask. (Default: %(default)s)")
    parser.add_argument("--iso_resolution", type=float, default=1.0,
                        help="Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
    parser.add_argument("--femur_side", type=int, default=1,
                        help="Set whether the left or right femur is to be analyzed. 1 = Left; 2 = Right. (default: %(default)s)")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Set output path for the N88 model file. (default: same as input image)")
    parser.add_argument("--elastic_exponent", type=float, default=2.29,
                        help="Sets the exponent (b) for power law: E=A(den)^b, used in converting material IDS/density to Elastic Moduli value. (default: %(default)s)")
    parser.add_argument("--elastic_Emax", type=float, default=10500,
                        help="Sets the coefficient (A) Elastic Modulus value for the power law: E=A(den)^b. (default: %(default)s [MPa])")
    parser.add_argument("--poissons_ratio", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_E", type=float, default=2500,
                        help="Sets the Elastic Modulus for PMMA caps in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_v", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the PMMA material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_thick", type=float, default=2.0,
                        help="Sets the minimum thickness for PMMA caps in the FE model. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=5000,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=-1.0,
                        help="Sets the applied displacement in [mm] to the FE model. (default: %(default)s [mm])")
    parser.add_argument("--reference_path", type=str, required=False, default=None,
                        help="Path to the reference vtk file for ICP registration. (default: None)")
    parser.add_argument("--pmma_yield_compression", type=float, default=None,
                        help="Sets the yield strength in compression for PMMA material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_yield_tension", type=float, default=None,
                        help="Sets the yield strength in tension for PMMA material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--bone_yield_compression", type=float, default=None,
                        help="Sets the yield strength in compression for bone material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--bone_yield_tension", type=float, default=None,
                        help="Sets the yield strength in tension for bone material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pad_vox", type=int, default=20,
                        help="Number of voxels to pad on all sides before femur pre-rotation and ICP transform. (default: %(default)s)")

    # Parse and display
    args = parser.parse_args()

    # Set default reference paths
    if args.reference_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, "../../dat")
        args.left_femur_reference = os.path.join(data_dir, "LT_FEMUR_SIDEWAYS_FALL_REF.vtk")
        args.right_femur_reference = os.path.join(data_dir, "RT_FEMUR_SIDEWAYS_FALL_REF.vtk")
    else:
        args.left_femur_reference = os.path.join(args.reference_path, "LT_FEMUR_SIDEWAYS_FALL_REF.vtk")
        args.right_femur_reference = os.path.join(args.reference_path, "RT_FEMUR_SIDEWAYS_FALL_REF.vtk")

    basename = remove_extension(os.path.basename(args.calibrated_image))

    if args.output_path is None:
        output_dir = os.path.dirname(args.calibrated_image)
    else:
        output_dir = args.output_path

    output_dir = os.path.abspath(output_dir)
    args.output_file = os.path.join(output_dir, f"{basename}.n88model")

    print(echo_arguments('ogoSidewaysFallFe', vars(args)))

    # Run program
    sidewaysFallFe(args)


if __name__ == '__main__':
    main()