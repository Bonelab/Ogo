#####
# ogo_l4_compression_fe.py
#
# This script sets up the L4 vertebral compression FE model from the density (K2HPO4)
# calibrated image. This script sets up the model for a L4 vertebra (including arch and
# pedicles). The analysis resamples the image to isotropic voxels, transforms
# the image, applies the bone mask and bins the data. It then creates the FE model for
# solving using FAIM (>v8.0, Numerics Solutions Ltd, Calgary, Canada - Steven  Boyd).
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# July 11, 2019
# Modifed to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import ogo.cli.calib.ogo_helper as ogo
import os
import sys
import argparse
import time
from datetime import date
from collections import OrderedDict
import vtk
import vtkbone


from ogo.util.echo_arguments import echo_arguments

##
# Reference File Associated with this script
l4_reference = "L4_BODY_SPINE_COMPRESSION_REF.vtk"


##
# Start script
def l4Compression_fe(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask

    mask_threshold = args.mask_threshold
    isotropic_voxel = args.isotropic_voxel
    elastic_exponent = args.elastic_exponent
    elastic_Emax = args.elastic_Emax
    poissons_ratio = args.poissons_ratio
    pmma_E = args.pmma_E
    pmma_v = args.pmma_v
    pmma_thick = args.pmma_thick
    pmma_mat_id = args.pmma_mat_id
    fe_displacement = args.fe_displacement


    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)
    script_name = sys.argv[0]
    N88_fileName = mask_basename.replace("_PERI_CORR.nii", "_L4_FE.n88model")

    ##
    # Message the input parameters to the terminal
    ogo.message("Image Path: %s" % image_pathname)
    ogo.message("Image File: %s" % image_basename)
    ogo.message("Mask Path: %s" % mask_pathname)
    ogo.message("Mask File: %s" % mask_basename)
    ogo.message("L4 Reference Model: %s" % l4_reference)
    ogo.message("Isotropic Voxel Size: %8.4f" % isotropic_voxel)
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
    ogo.message("Resampling the Input image and bone mask to isotropic...")
    image_resample = ogo.imageResample(imageData, isotropic_voxel)
    mask_resample = ogo.imageResample(maskThres, isotropic_voxel)

    ##
    # Extracting vertebral body from full vertebra
    ogo.message("Extracting vertebral body from full vertebra...")
    vertebral_body_image, vertebral_body_mask = ogo.vertebralBodyExtract(image_resample, mask_resample)

    ##
    # Align the input femur with the reference model
    ogo.message("Aligning input with reference model...")
    mask_surface = ogo.marchingCubes(vertebral_body_mask)
    ref_poly = ogo.readPolyData(l4_reference)

    icp = ogo.iterativeClosestPoint(ref_poly, mask_surface)

    ogo.message("Applying the transformation to the image and mask...")
    image_trans = ogo.applyTransform(vertebral_body_image, icp)
    mask_trans = ogo.applyTransform(vertebral_body_mask, icp)

    # ogo.message("Writing out temp images...")
    # ogo.writeNii(image_trans, "temp_image.nii", image_pathname)
    # ogo.writeNii(mask_trans, "temp_mask.nii", image_pathname)
    #
    # l4_reference_image = "/Users/andrew.michalski/Development/Ogo/L4_BODY_SPINE_COMPRESSION_REF.nii"
    # ogo.message("Applying SimpleITK registration of the images...")
    # if femur_side == 1:
    #     ogo.finalRegistration(left_femur_reference_image)
    # elif femur_side == 2:
    #     ogo.finalRegistration(right_femur_reference_image)
    # ogo.message("Reading temp images...")
    # image_trans = ogo.readNii("temp_image2.nii")
    # mask_trans = ogo.readNii("temp_mask2.nii")
    # os.remove("temp_image2.nii")
    # os.remove("temp_mask2.nii")

    # Replace any negative values less than -31 to be equivalent to -31.
    # -31 is used as that converts to a minumum elastic modulus value of 0.1 MPa
    # K2HPO4 den = -31 mg/cc => Ash den = 6 mg/cc => E = 0.1 MPa
    image_thres = ogo.bmd_preprocess(image_trans, -31)
    image_ash = ogo.bmd_K2hpo4ToAsh(image_thres)

    ##
    # Set up the FE model.
    ogo.message("Setting up the Finite Element Model...")
    # Cast the image to Short to "round" float values to nearest whole number
    cast_image = ogo.cast2short(image_ash)
    cast_mask = ogo.cast2unsignchar(mask_trans)

    # Apply the mask to the bone
    ogo.message("Applying the bone mask to the image...")
    bone_image = ogo.applyMask(cast_image, cast_mask)

    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn = ogo.imageConnectivity(bone_image)

    change = ogo.changeInfo(conn)

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh = ogo.Image2Mesh(change)

    ##
    # Set up the Material Table
    ogo.message("Setting up the Finite Element Material Table...")
    material_table = ogo.materialTable(
        mesh,
        poissons_ratio,
        elastic_Emax,
        elastic_exponent,
        pmma_mat_id,
        pmma_E,
        pmma_v
    )

    ##
    # Create the preliminary FE model
    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()
    model_bounds = model.GetBounds()
    ogo.message("Model Bounds: %s" % str(model_bounds))

    # Define support vectors for later application
    z_center = (model_bounds[4] + model_bounds[5]) / 2
    top_support_vector = (0, 0, 1)
    bottom_support_vector = (0, 0, -1)

    ogo.message("Determining the superior vertebral body...")
    superior_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[2]),
        float(model_bounds[3]),
        float(model_bounds[5] - model_bounds[5]/25),
        float(model_bounds[5])
    )

    ogo.message("Extracted Superior Bounds: %s" % str(superior_bounds))
    superior_model = ogo.extractBox(superior_bounds, model)
    superior_model_bounds = superior_model.GetBounds()
    ogo.message("Extracted Superior Model Bounds: %s" % str(superior_model_bounds))

    superior_model_bounds = (
        superior_model_bounds[0],
        superior_model_bounds[1] - 1,
        superior_model_bounds[2],
        superior_model_bounds[3] - 1,
        superior_model_bounds[4],
        superior_model_bounds[5] - 1
    )

    superior_element_dims = (
        superior_model_bounds[1] - superior_model_bounds[0],
        superior_model_bounds[3] - superior_model_bounds[2],
        superior_model_bounds[5] - superior_model_bounds[4]
    )

    ##
    # Extract inferior vertebral body
    ogo.message("Determining the inferior vertebral body...")
    inferior_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[2]),
        float(model_bounds[3]),
        float(model_bounds[4]),
        float(model_bounds[4] + model_bounds[4]/25)
    )
    ogo.message("Extracted Inferior Bounds: %s" % str(inferior_bounds))
    inferior_model = ogo.extractBox(inferior_bounds, model)
    inferior_model_bounds = inferior_model.GetBounds()
    ogo.message("Extracted Inferior Model Bounds: %s" % str(inferior_model_bounds))

    inferior_model_bounds = (
        inferior_model_bounds[0],
        inferior_model_bounds[1] - 1,
        inferior_model_bounds[2],
        inferior_model_bounds[3] - 1,
        inferior_model_bounds[4],
        inferior_model_bounds[5] - 1
    )

    inferior_element_dims = (
        inferior_model_bounds[1] - inferior_model_bounds[0],
        inferior_model_bounds[3] - inferior_model_bounds[2],
        inferior_model_bounds[5] - inferior_model_bounds[4]
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
    # Creates the femoral head pmma cap
    ogo.message("Creating Superior PMMA Cap...")
    superiorPMMA = ogo.superiorVertebralPMMA(
        superior_model_bounds,
        spacing,
        origin,
        inval,
        outval,
        thickness,
        pmma_mat_id
    )

    ##
    # Creates the femoral head pmma cap
    ogo.message("Creating Inferior PMMA Cap...")
    inferiorPMMA = ogo.inferiorVertebralPMMA(
        inferior_model_bounds,
        spacing,
        origin,
        inval,
        outval,
        thickness,
        pmma_mat_id
    )

    ##
    # combines the pmma cap images with the model image.
    ogo.message("Combine PMMA Cap Images with Model Image...")
    combinedImage = ogo.combineImageData_VC(bone_image, superiorPMMA, inferiorPMMA, pmma_mat_id)

    ##
    # Mesh the final image and create Finite Element Model
    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn2 = ogo.imageConnectivity(combinedImage)

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh2 = ogo.Image2Mesh(conn2)

    ##
    # Create the final FE model
    ogo.message("Constructing the Finite Element Model...")
    model2 = ogo.applyTestBase(mesh2, material_table)
    model2.ComputeBounds()
    model2_bounds = model2.GetBounds()
    ogo.message("Model 2 Bounds: %s" % str(model2_bounds))
    z_center2 = (model2_bounds[4] + model2_bounds[5]) / 2

    ##
    # Determine the superior PMMA cap nodes
    ogo.message("Determining the Superior PMMA Cap nodes...")
    sv_pmma_bounds = (
        model2_bounds[0],
        model2_bounds[1],
        model2_bounds[2],
        model2_bounds[3],
        model2_bounds[5]-1,
        model_bounds[5]
    )
    ogo.message("Superior PMMA Cap Bounds: %s" % str(sv_pmma_bounds))
    sv_pmma_model = ogo.extractBox(sv_pmma_bounds, model2)
    sv_pmma_model_bounds = sv_pmma_model.GetBounds()
    ogo.message("Superior PMMA Cap Model Bounds: %s" % str(sv_pmma_model_bounds))

    # Get the visible nodes on superiot pmma cap
    sv_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        sv_pmma_visible_node_IDS,
        sv_pmma_model,
        top_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Superior PMMA Cap."
        % sv_pmma_visible_node_IDS.GetNumberOfTuples())
    sv_pmma_visible_node_IDS.SetName("Superior_PMMA_Nodes")
    model2.AddNodeSet(sv_pmma_visible_node_IDS)

    ##
    # Determine the inferior PMMA cap nodes
    ogo.message("Determining the Inferior PMMA Cap nodes...")
    iv_pmma_bounds = (
        model2_bounds[0],
        model2_bounds[1],
        model2_bounds[2],
        model2_bounds[3],
        model2_bounds[4],
        model_bounds[4]+1
    )
    ogo.message("Inferior PMMA Cap Bounds: %s" % str(iv_pmma_bounds))
    iv_pmma_model = ogo.extractBox(iv_pmma_bounds, model2)
    iv_pmma_model_bounds = iv_pmma_model.GetBounds()
    ogo.message("Inferior PMMA Cap Model Bounds: %s" % str(iv_pmma_model_bounds))

    # Get the visible nodes on superiot pmma cap
    iv_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        iv_pmma_visible_node_IDS,
        iv_pmma_model,
        bottom_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Inferior PMMA Cap."
        % iv_pmma_visible_node_IDS.GetNumberOfTuples())
    iv_pmma_visible_node_IDS.SetName("Inferior_PMMA_Nodes")
    model2.AddNodeSet(iv_pmma_visible_node_IDS)

    ##
    # Apply Boundary Conditions to PMMA caps
    # - Superior PMMA Cap BCs
    ogo.message("Applying boundary conditions to the Superior PMMA cap...")
    model2.ApplyBoundaryCondition(
        "Superior_PMMA_Nodes",
        vtkbone.vtkboneConstraint.SENSE_Z,
        fe_displacement,
        "Superior_PMMA_Z_Displacement"
    )

    # - Inferior PMMA Cap BCs
    ogo.message("Applying boundary conditions to the Inferior PMMA cap...")
    model2.ApplyBoundaryCondition(
        "Inferior_PMMA_Nodes",
        vtkbone.vtkboneConstraint.SENSE_Z,
        0,
        "Inferior_PMMA_Z_Displacement"
    )

    ##
    # Setting up the post-processing Parameters
    ogo.message("Setting up the Post Processing Parameters...")
    info = model2.GetInformation()
    pp_node_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_NODE_SETS()
    pp_elem_sets_key = vtkbone.vtkboneSolverParameters.POST_PROCESSING_ELEMENT_SETS()

    for setname in ["Superior_PMMA_Nodes", "Inferior_PMMA_Nodes"]:
        # Specify that these sets should be used for post-processing
        pp_node_sets_key.Append(info, setname)
        elementSet = model2.GetAssociatedElementsFromNodeSet(setname)
        model2.AddElementSet(elementSet)
        pp_elem_sets_key.Append(info, setname)

    model2.AppendHistory(
        "Created by %s version %s." % (script_name, script_version))

    ##
    # Write out n88model file
    ogo.message("Writing out n88model file: %s" % N88_fileName)
    ogo.writeN88Model(model2, N88_fileName, image_pathname)

    ##
    ogo.message("Done writing n88model.")

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()

def main():
    description = '''
    This script sets up the L4 vertebral compression FE model from the
        density (K2HPO4) calibrated image. This script sets up the model for a L4 vertebra
        (including arch and pedicles). The analysis resamples the image to isotropic voxels,
        transforms the image, applies the bone mask and bins the data. It then creates the FE
        model for solving using FAIM (v8.1, Numerics Solutions Ltd, Calgary, Canada - Steven
        Boyd).

        Input: Calibrated K2HPO4 Image (*.nii), Bone Mask (*_MASK.nii)

        Optional Parameters:
        1) Mask Threshold
        2) Isotropic resample voxel size
        3) Power-law exponent
        4) Power-law coefficient
        5) Bone Poissons ratio
        6) PMMA Elastic Modulus
        7) PMMA Poissons ratio
        8) PMMA thickness
        9) PMMA material ID
        10) FE displacement

        Output: N88 Model (*.n88model)
    '''


    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="ogoL4CompressionFe",
        description=description
    )

    parser.add_argument("calibrated_image",
        help = "*_K2HPO4.nii image file")
    parser.add_argument("bone_mask",
        help = "*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type = int,
                        default = 7,
                        help = "Set the threshold value to extract the bone of interest from the mask. (Default: %(default)s)")
    parser.add_argument("--isotropic_voxel", type = float, default = 1.0,
                        help = "Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
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
    parser.add_argument("--pmma_thick", type=float, default=3.0,
                        help="Sets the minimum thickness for PMMA caps in the FE model. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=5000,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=-1.0,
                        help="Sets the applied displacement in [mm] to the FE model. (default: %(default)s [mm])")

    # Parse and display
    args = parser.parse_args()
    
    print(echo_arguments('l4_compression_fe', vars(args)))

    # Run program
    l4Compression_fe(args)

if __name__ == '__main__':
    main()