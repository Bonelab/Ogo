#####
# ogo_sideways_fall_fe.py
#
# This script sets up the sideways fall FE model on the hip from the density (K2HPO4)
# calibrated image. This script sets up the model for either a left or right femur, as
# specified by the user. The analysis resamples the image to isotropic voxels, transforms
# the image, applies the bone mask and bins the data. It then creates the FE model for
# solving using FAIM (>v8.0, Numerics Solutions Ltd, Calgary, Canada - Steven  Boyd).
#
#####
#
# Andrew Michalski
# University of Calgary
# Biomedical Engineering Graduate Program
# April 29, 2019
# Modified to Py3: March 25, 2020
#####

script_version = 1.0

##
# Import the required modules
import ogo.util.Helper as ogo
import math
import os
import sys
import argparse
import time
from datetime import date
from collections import OrderedDict
import vtk
import vtkbone

from ogo.fea.boundary import (
    generate_fixture_cap_vtk,
    label_foreground_in_bounds_vtk,
    should_smooth_resampled_mask,
    smooth_binary_mask_vtk,
)
from ogo.fea.femur import (
    DEFAULT_CORTICAL_LABEL,
    DEFAULT_FEMUR_FE_DISPLACEMENT,
    DEFAULT_FEMUR_CUT_MODE,
    DEFAULT_FEMUR_INPUT_MARGIN_MM,
    DEFAULT_FEMUR_ISO_RESOLUTION_MM,
    DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
    DEFAULT_FEMUR_SHAFT_LENGTH_MM,
    DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
    DEFAULT_PMMA_INTRUSION_MM,
    DEFAULT_PMMA_THICKNESS_MM,
    DEFAULT_TRABECULAR_LABEL,
    DISTAL_FEMUR_NODE_SET,
    FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM,
    FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM,
    FEMORAL_HEAD_NODE_SET,
    GREATER_TROCHANTER_NODE_SET,
    SIDEWAYS_FALL_NODE_SETS,
    cortical_compartment_mask,
    expand_xz_footprint,
    femur_z_coverage,
    flat_crop_vtk_image_below_z,
    mirror_polydata_x,
    pad_vtk_images_to_foreground_margin,
    side_rotation,
    sideways_fall_output_name,
    standardize_femur_shaft_length,
    swap_xz_footprint,
)
from ogo.fea.materials import build_femur_material_table
from ogo.fea.model import append_postprocessing_sets, find_nodes_on_coordinate_plane, write_model
from ogo.util.echo_arguments import echo_arguments


def remove_extension(filename):
    while True:
        filename, ext = os.path.splitext(filename)
        if not ext:
            break
    return filename


##
# Start script
def sidewaysFallFe(args):
    ogo.message("Start of Script...")

    ##
    # Collect the input arguments
    image = args.calibrated_image
    mask = args.bone_mask
    compartment_mask = args.compartment_mask

    mask_threshold = args.mask_threshold
    iso_resolution = args.iso_resolution
    femur_side = args.femur_side
    poissons_ratio = args.poissons_ratio
    pmma_E = args.pmma_E
    pmma_v = args.pmma_v
    pmma_thick = args.pmma_thick
    pmma_intrusion = args.pmma_intrusion
    pmma_mat_id = args.pmma_mat_id
    fe_displacement = args.fe_displacement
    left_femur_reference = args.left_femur_reference
    right_femur_reference = args.right_femur_reference
    output_file = args.output_file
    pmma_yield_compression = args.pmma_yield_compression
    pmma_yield_tension = args.pmma_yield_tension
    femur_shaft_length = args.femur_shaft_length
    femur_cut_mode = args.femur_cut_mode
    femur_lesser_trochanter_distal_offset = args.femur_lesser_trochanter_distal_offset
    femur_lesser_trochanter_distal_offset_percent = args.femur_lesser_trochanter_distal_offset_percent
    femur_input_margin = args.femur_input_margin
    cortical_label = args.cortical_label
    trabecular_label = args.trabecular_label
    mask_smoothing_spacing_threshold = args.mask_smoothing_spacing_threshold



    ##
    # Determine image locations and names of files
    image_pathname = os.path.dirname(image)
    image_basename = os.path.basename(image)
    mask_pathname = os.path.dirname(mask)
    mask_basename = os.path.basename(mask)
    script_name = sys.argv[0]
    
    
    try:
        N88_fileName = sideways_fall_output_name(output_file, femur_side)
        rot_z = side_rotation(femur_side)
    except ValueError:
        ogo.message("Femur side not recognized. Terminating...")
        sys.exit()


    ##
    # Message the input parameters to the terminal
    ogo.message("Image Path: %s" % image_pathname)
    ogo.message("Image File: %s" % image_basename)
    ogo.message("Mask Path: %s" % mask_pathname)
    ogo.message("Mask File: %s" % mask_basename)
    ogo.message("Isotropic Voxel Size: %8.4f" % iso_resolution)
    if femur_side == 1:
        ogo.message("Femur Side for Model: Left")
    elif femur_side == 2:
        ogo.message("Femur Side for Model: Right")
    else:
        ogo.message("Femur side not recognized. Terminating...")
        sys.exit()
    ogo.message("Trabecular Elastic Function: %s" % (args.elastic_E_func or "default_E"))
    ogo.message("Cortical Elastic Function: %s" % (args.cort_elastic_E_func or args.elastic_E_func or "default_E"))
    ogo.message("Bone Poissons Ratio: %1.1f" % poissons_ratio)
    ogo.message("PMMA Elastic Modulus [MPa]: %8.4f" % pmma_E)
    ogo.message("PMMA Poissons Ration: %1.1f" % pmma_v)
    ogo.message("PMMA Thickness [mm]: %8.4f" % pmma_thick)
    ogo.message("PMMA Intrusion [mm]: %8.4f" % pmma_intrusion)
    ogo.message("PMMA Material ID: %d" % pmma_mat_id)
    ogo.message("Applied Displacement [mm]: %8.4f" % fe_displacement)
    ogo.message("Femur Cut Mode: %s" % femur_cut_mode)
    if compartment_mask is not None:
        ogo.message("Compartment Mask: %s" % compartment_mask)
        ogo.message("Cortical Label: %d" % cortical_label)
        ogo.message("Trabecular Label: %d" % trabecular_label)
    if femur_cut_mode == "fixed_length":
        ogo.message("Retained Proximal Femur Length [mm]: %8.4f" % femur_shaft_length)
    elif femur_lesser_trochanter_distal_offset_percent is not None:
        ogo.message(
            "Lesser Trochanter Distal Cut Offset [%% greater-lesser distance]: %8.4f"
            % femur_lesser_trochanter_distal_offset_percent
        )
    else:
        ogo.message(
            "Lesser Trochanter Distal Cut Offset [mm]: %8.4f"
            % femur_lesser_trochanter_distal_offset
        )
    ogo.message("Input Foreground Safety Margin [mm]: %8.4f" % femur_input_margin)

    ##
    # Read input image
    ogo.message("Reading calibrated image...")
    imageData = ogo.readNii(image)
    input_spacing = imageData.GetSpacing()

    ##
    # Read bone mask
    ogo.message("Reading bone mask...")
    maskData = ogo.readNii(mask)
    maskThres = ogo.maskThreshold(maskData, mask_threshold)
    compartmentData = None
    if compartment_mask is not None:
        ogo.message("Reading trabecular/cortical compartment mask...")
        compartmentData = ogo.readNii(compartment_mask)

    images_to_pad = [imageData, maskThres] + ([compartmentData] if compartmentData is not None else [])
    padded_images, padding = pad_vtk_images_to_foreground_margin(
        images_to_pad,
        maskThres,
        margin_mm=femur_input_margin,
        constants=[0] * len(images_to_pad),
    )
    imageData = padded_images[0]
    maskThres = padded_images[1]
    if compartmentData is not None:
        compartmentData = padded_images[2]
    if any(padding["lower"]) or any(padding["upper"]):
        ogo.message(
            "Padded input image extent by lower=%s upper=%s voxels for FE transform safety."
            % (padding["lower"], padding["upper"])
        )
    else:
        ogo.message("Input image already has sufficient foreground safety margin.")

    ##
    # Pre-rotate the image and mask for better alignment
    image_rot, mask_rot = ogo.preRotateImage(imageData, maskThres, rot_z)
    compartment_rot = None
    if compartmentData is not None:
        _, compartment_rot = ogo.preRotateImage(imageData, compartmentData, rot_z)

    ##
    # Align the input femur with the reference model
    ogo.message("Aligning input with reference model...")
    mask_surface = ogo.marchingCubes(mask_rot)
    if femur_side == 1:
        ref_poly = ogo.readPolyData(left_femur_reference)
    elif femur_side == 2:
        if os.path.exists(right_femur_reference):
            ref_poly = ogo.readPolyData(right_femur_reference)
        else:
            ogo.message(
                "Right femur reference not found; mirroring left reference in x:",
                right_femur_reference,
            )
            ref_poly = mirror_polydata_x(ogo.readPolyData(left_femur_reference))
    else:
        print("Error: Femur Side not defined. Terminating...")
        sys.exit()

    icp = ogo.iterativeClosestPoint(ref_poly, mask_surface)
    # ogo.message("Transform:\n %s" % str(icp))

    ogo.message("Applying the transformation and isotropic resampling to the image and mask...")
    image_trans = ogo.transformResample(image_rot, icp, iso_resolution)
    mask_trans = ogo.labelTransformResample(mask_rot, icp, iso_resolution)
    compartment_trans = (
        ogo.labelTransformResample(compartment_rot, icp, iso_resolution)
        if compartment_rot is not None
        else None
    )
    smooth_resampled_masks = should_smooth_resampled_mask(
        input_spacing,
        mask_smoothing_spacing_threshold,
    )
    if smooth_resampled_masks:
        ogo.message(
            "smoothing resampled femur mask because input spacing "
            f"{input_spacing} exceeds {mask_smoothing_spacing_threshold} mm..."
        )
        mask_trans = smooth_binary_mask_vtk(mask_trans, close_iter=1, open_iter=1)
    else:
        ogo.message(
            "skipping femur mask smoothing because input spacing "
            f"{input_spacing} is <= {mask_smoothing_spacing_threshold} mm in all dimensions..."
        )

    ogo.message("Applying flat distal femur crop in reference coordinates...")
    try:
        image_trans, mask_trans, shaft_crop = standardize_femur_shaft_length(
            image_trans,
            mask_trans,
            reference_bounds=ref_poly.GetBounds(),
            retained_length_mm=femur_shaft_length,
            cut_mode=femur_cut_mode,
            lesser_trochanter_distal_offset_mm=femur_lesser_trochanter_distal_offset,
            lesser_trochanter_distal_offset_percent=femur_lesser_trochanter_distal_offset_percent,
        )
    except ValueError as exc:
        ogo.message(str(exc))
        sys.exit(1)
    retained_length_mm = shaft_crop["retained_length_mm"]
    ogo.message(
        "Distal cut z=%8.4f; input mask z coverage [%8.4f, %8.4f]; retained length=%8.4f"
        % (
            shaft_crop["cut_z"],
            shaft_crop["mask_z_min"],
            shaft_crop["mask_z_max"],
            retained_length_mm,
        )
    )
    for warning in shaft_crop.get("warnings", []):
        ogo.message("WARNING: Femur crop QC: %s" % warning)
    if femur_cut_mode == "lesser_trochanter":
        ogo.message(
            "Detected greater trochanter z=%8.4f; lesser trochanter z=%8.4f"
            % (shaft_crop["greater_trochanter_z"], shaft_crop["lesser_trochanter_z"])
        )
    if compartment_trans is not None:
        compartment_trans = flat_crop_vtk_image_below_z(compartment_trans, shaft_crop["cut_z"])

    # ogo.message("Writing out temp images...")
    # ogo.writeNii(image_trans, "temp_image.nii", image_pathname)
    # ogo.writeNii(mask_trans, "temp_mask.nii", image_pathname)
    #
    # left_femur_reference_image = "/Users/andrew.michalski/Development/Ogo/LT_FEMUR_SIDEWAYS_FALL_REF.nii"
    # right_femur_reference_image = "/Users/andrew.michalski/Development/Ogo/RT_FEMUR_SIDEWAYS_FALL_REF.nii"
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

    # replace any negative values less than -31 to be equivalent to -31.
    # -31 is used as that converts to a minumum elastic modulus value of 0.1 MPa
    # K2HPO4 den = -31 mg/cc => Ash den = 6 mg/cc => E = 0.1 MPa
    image_thres = ogo.bmd_preprocess(image_trans, -31)
    image_ash = ogo.bmd_K2hpo4ToAsh(image_thres)

    cortical_mask = None
    if compartment_trans is not None:
        cortical_mask = cortical_compartment_mask(
            compartment_trans,
            cortical_label=cortical_label,
            trabecular_label=trabecular_label,
        )
        if smooth_resampled_masks:
            ogo.message("smoothing derived cortical compartment mask...")
            cortical_mask = smooth_binary_mask_vtk(cortical_mask, close_iter=1, open_iter=1)
    binned_image, bin_centers = ogo.density2materialID(
        image_ash,
        n_bins=128,
        cort_mask=cortical_mask,
    )

    ##
    # Set up the FE model.
    ogo.message("Setting up the Finite Element Model...")
    # Cast the image to Short to "round" float values to nearest whole number
    cast_image = ogo.cast2short(binned_image)
    cast_mask = ogo.cast2unsignchar(mask_trans)


    # Apply the mask to the bone
    ogo.message("Applying the bone mask to the image...")
    bone_image = ogo.applyMask(cast_image, cast_mask)

    # Ensure connectivity
    ogo.message("Performing Image connectivity...")
    conn = ogo.imageConnectivity(bone_image)


    change = ogo.changeInfo(conn)
    model_z_min, model_z_max = femur_z_coverage(change)
    model_cut_z = model_z_max - retained_length_mm
    if model_z_min > model_cut_z:
        if model_z_min - model_cut_z <= 2.0:
            model_cut_z = model_z_min
        else:
            ogo.message(
                "Femur scan does not include enough distal shaft after model-grid alignment "
                "(required retained length %.1f mm, available approximately %.1f mm)."
                % (retained_length_mm, max(0.0, model_z_max - model_z_min))
            )
            sys.exit(1)
    change = flat_crop_vtk_image_below_z(change, model_cut_z)
    ogo.message(
        "Applied model-grid distal shaft cut z=%8.4f; model mask z coverage [%8.4f, %8.4f]"
        % (model_cut_z, model_z_min, model_z_max)
    )

    # Convert image data to hexahedral elements
    ogo.message("Meshing image data to elements...")
    mesh = ogo.Image2Mesh(change)

    ##
    # Set up the Material Table
    ogo.message("Setting up the Finite Element Material Table...")
    material_table = build_femur_material_table(
        bin_centers,
        n_bins=128,
        elastic_E_func=args.elastic_E_func,
        yield_comp_func=args.yield_comp_func,
        yield_tens_func=args.yield_tens_func,
        cort_elastic_E_func=args.cort_elastic_E_func,
        cort_yield_comp_func=args.cort_yield_comp_func,
        cort_yield_tens_func=args.cort_yield_tens_func,
        cort_poissons_ratio=args.cort_poissons_ratio,
        poissons_ratio=poissons_ratio,
        pmma_mat_id=pmma_mat_id,
        pmma_E=pmma_E,
        pmma_v=pmma_v,
        pmma_yield_tension=pmma_yield_tension,
        pmma_yield_compression=pmma_yield_compression,
        include_cortical=cortical_mask is not None,
    )
    ##
    # Create the preliminary FE model
    ogo.message("Constructing the Finite Element Model...")
    model = ogo.applyTestBase(mesh, material_table)
    model.ComputeBounds()
    model_bounds = model.GetBounds()
    ogo.message("Model Bounds: %s" % str(model_bounds))

    # Define support vectors for later application
    z_center = (model_bounds[4] + model_bounds[5] / 2)
    top_support_vector = (0, 1, 0)
    bottom_support_vector = (0, -1, 0)
    side_support_vector = (0, 0, -1)

    # The shaft is cut on a flat z plane, so use that cut-plane normal for the distal BC.
    side_support_vector = (0, 0, -1)

    ogo.message("Determining Femoral Head ROI...")
    femoral_head_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[2]),
        float(model_bounds[2] + (model_bounds[3]-model_bounds[2]) / 5),
        float(model_bounds[4]),
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

    femoral_head_element_dims = (
        femoral_head_model_bounds[1] - femoral_head_model_bounds[0],
        femoral_head_model_bounds[3] - femoral_head_model_bounds[2],
        femoral_head_model_bounds[5] - femoral_head_model_bounds[4]
    )

    ##
    # Extract Greater Trochater ROI
    ogo.message("Determining Greater Trochater ROI...")
    greater_trochanter_bounds = (
        float(model_bounds[0]),
        float(model_bounds[1]),
        float(model_bounds[3]-(model_bounds[3]-model_bounds[2])/20),
        float(model_bounds[3]),
        float(model_bounds[4]),
        float(model_bounds[5])
    )
    ogo.message("Extracted Greater Trochanter Bounds: %s" % str(greater_trochanter_bounds))
    greater_trochanter_model = ogo.extractBox(greater_trochanter_bounds, model)
    greater_trochanter_model_bounds = greater_trochanter_model.GetBounds()
    ogo.message("Extracted Greater Trochanter Model Bounds: %s"
        % str(greater_trochanter_model_bounds))

    greater_trochanter_model_bounds = (
        greater_trochanter_model_bounds[0],
        greater_trochanter_model_bounds[1] - 1,
        greater_trochanter_model_bounds[2],
        greater_trochanter_model_bounds[3] - 1,
        greater_trochanter_model_bounds[4],
        greater_trochanter_model_bounds[5] - 1
    )

    greater_trochanter_element_dims = (
        greater_trochanter_model_bounds[1] - greater_trochanter_model_bounds[0],
        greater_trochanter_model_bounds[3] - greater_trochanter_model_bounds[2], greater_trochanter_model_bounds[5] - greater_trochanter_model_bounds[4]
    )

    # Thickness to be added to PMMA cap.
    thickness = max(1, int(math.ceil(pmma_thick / change.GetSpacing()[1])))
    intrusion = max(0, int(math.ceil(pmma_intrusion / change.GetSpacing()[1])))

    ##
    # Creates the femoral head pmma cap
    ogo.message("Creating Femoral Head PMMA Cap from contact ROI...")
    femoral_head_fixture_bounds = expand_xz_footprint(
        swap_xz_footprint(femoral_head_model_bounds),
        x_extension_mm=FEMORAL_HEAD_FIXTURE_WIDTH_EXTENSION_MM,
        z_extension_mm=FEMORAL_HEAD_FIXTURE_LONG_AXIS_EXTENSION_MM,
    )
    ogo.message("Femoral Head PMMA Fixture Bounds: %s" % str(femoral_head_fixture_bounds))
    femoral_head_contact = label_foreground_in_bounds_vtk(change, femoral_head_fixture_bounds, label_value=1)
    femoralHeadPMMA = generate_fixture_cap_vtk(
        femoral_head_contact,
        label_value=1,
        axis="y",
        direction="down",
        thickness=thickness,
        shape="box",
        intrusion=intrusion,
        crop_to_contact=True,
        output_value=pmma_mat_id,
    )

    ##
    # Creates the greater trochanter pmma cap
    ogo.message("Creating Greater Trochanter PMMA Cap from contact ROI...")
    greater_trochanter_contact = label_foreground_in_bounds_vtk(change, greater_trochanter_model_bounds, label_value=1)
    greaterTrochanterPMMA = generate_fixture_cap_vtk(
        greater_trochanter_contact,
        label_value=1,
        axis="y",
        direction="up",
        thickness=thickness,
        shape="round",
        intrusion=intrusion,
        output_value=pmma_mat_id,
    )

    ##
    # combines the pmma cap images with the model image.
    ogo.message("Combine PMMA Cap Images with Model Image...")
    combinedImage = ogo.combineImageData_SF(change, femoralHeadPMMA, greaterTrochanterPMMA, pmma_mat_id)

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
    z_center2 = (model2_bounds[4] + model2_bounds[5] / 2)

    ##
    # Determine Femoral Head PMMA Cap support nodes.
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

    # Get the visible nodes of the femoral head pmma cap
    fh_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        fh_pmma_visible_node_IDS,
        fh_pmma_model,
        bottom_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Femoral Head PMMA Cap."
        % fh_pmma_visible_node_IDS.GetNumberOfTuples())
    fh_pmma_visible_node_IDS.SetName(FEMORAL_HEAD_NODE_SET)
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

    # Get the visible nodes on the Greater Trochanter PMMA Cap
    gt_pmma_visible_node_IDS = vtk.vtkIdTypeArray()
    vtkbone.vtkboneNodeSetsByGeometry.FindNodesOnVisibleSurface(
        gt_pmma_visible_node_IDS,
        gt_pmma_model,
        top_support_vector,
        -1
    )

    ogo.message("-- found %d visible exterior nodes on Greater Trochanter PMMA Cap."
        % gt_pmma_visible_node_IDS.GetNumberOfTuples())

    gt_pmma_visible_node_IDS.SetName(GREATER_TROCHANTER_NODE_SET)
    model2.AddNodeSet(gt_pmma_visible_node_IDS)

    ##
    # Distal Femur (df) support nodes. The shaft is cut on a flat z plane, so
    # the fixed support is the cut face itself rather than a thick distal band.
    ogo.message("Determining distal femur nodes...")
    distal_cut_z = model2_bounds[4]
    ogo.message("Distal Femur Cut Plane z: %8.4f" % distal_cut_z)
    df_visible_node_IDS = find_nodes_on_coordinate_plane(model2, "z", distal_cut_z)

    ogo.message("-- found %d visible exterior nodes on Distal Femur."
        % df_visible_node_IDS.GetNumberOfTuples())

    df_visible_node_IDS.SetName(DISTAL_FEMUR_NODE_SET)
    model2.AddNodeSet(df_visible_node_IDS)

    ##
    # Apply Boundary conditions to PMMA caps at specific sites
    ogo.message("Applying displacement boundary condition to Femoral Head PMMA cap...")
    model2.ApplyBoundaryCondition(
        FEMORAL_HEAD_NODE_SET,
        vtkbone.vtkboneConstraint.SENSE_Y,
        -fe_displacement,
        "top_displacement")

    ogo.message("Constraining Greater Trochanter PMMA cap in loading direction...")
    model2.ApplyBoundaryCondition(
        GREATER_TROCHANTER_NODE_SET,
        vtkbone.vtkboneConstraint.SENSE_Y,
        0,
        "bottom_fixed_y_PMMA")

    for SENSE, label in (
        (vtkbone.vtkboneConstraint.SENSE_X, "x"),
        (vtkbone.vtkboneConstraint.SENSE_Z, "z"),
    ):
        ogo.message("Constraining distal femur rigid-body motion...")
        model2.ApplyBoundaryCondition(
            DISTAL_FEMUR_NODE_SET,
            SENSE,
            0,
            f"bottom_fixed_{label}")

    ##
    # Post Processing parameters
    ogo.message("Setting up Post Processing Parameters...")
    append_postprocessing_sets(
        model2,
        SIDEWAYS_FALL_NODE_SETS,
    )

    model2.AppendHistory(
        "Created by %s version %s." % (script_name, script_version))

    ##
    # Write out n88model file
    ogo.message("Writing out n88model file: %s" % N88_fileName)
    write_model(model2, N88_fileName)

    ##
    ogo.message("Done writing n88model.")

    ##
    # End of script
    ogo.message("End of Script.")
    sys.exit()


def main():
     # Setup description
    description='''
This script sets up the sideways fall FE model on the hip from the
        density (K2HPO4) calibrated image. This script sets up the model for either a left or
        right femur, as specified by the user. The analysis resamples the image to isotropic
        voxels, transforms the image, applies the bone mask and bins the data. It then
        creates the FE model for solving using FAIM (v8.1, Numerics Solutions Ltd, Calgary,
        Canada - Steven  Boyd).

        Input: Calibrated K2HPO4 Image (*.nii), Bone Mask (*_MASK.nii)

        Optional Parameters:
        1) Isotropic resample voxel size
        2) Left or Right Femur Bone Mask

        Output: N88 Model (*.n88model)

'''
    

    # Setup argument parsing
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prog="OgoSidewaysFallFe",
        description=description
    )
    

    parser.add_argument("calibrated_image",
        help = "*_K2HPO4.nii image file")
    parser.add_argument("bone_mask",
        help = "*_MASK.nii mask image of bone")

    parser.add_argument("--mask_threshold", type = int,
                        default = 1,
                        help = "Set the threshold value to extract the bone of interest from the mask. (Default: %(default)s)")
    parser.add_argument("--iso_resolution", type = float, default = DEFAULT_FEMUR_ISO_RESOLUTION_MM,
                        help = "Set the isotropic voxel size [in mm]. (default: %(default)s [mm])")
    parser.add_argument("--mask_smoothing_spacing_threshold", type=float, default=DEFAULT_FEMUR_MASK_SMOOTHING_SPACING_THRESHOLD_MM,
                        help="Smooth resampled femur masks only when an input spacing dimension exceeds this value. (default: %(default)s [mm])")
    parser.add_argument("--femur_side", type = int, default = 1,
                        help = "Set whether the left or right femur is to be analyzed. 1 = Left; 2 = Right. (default: %(default)s)")
    parser.add_argument("--output_path", type=str, default=None,
                        help="Set output path for the N88 model file. (default: same as input image)")
    parser.add_argument("--poissons_ratio", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--elastic_E_func", type=str, default=None,
                        help="Function name for trabecular bone Young's modulus. (default: default_E)")
    parser.add_argument("--yield_comp_func", type=str, default=None,
                        help="Function name for trabecular compression yield. (default: none)")
    parser.add_argument("--yield_tens_func", type=str, default=None,
                        help="Function name for trabecular tension yield. (default: none)")
    parser.add_argument("--cort_elastic_E_func", type=str, default=None,
                        help="Function name for cortical bone Young's modulus. Defaults to --elastic_E_func.")
    parser.add_argument("--cort_yield_comp_func", type=str, default=None,
                        help="Function name for cortical compression yield. Defaults to --yield_comp_func.")
    parser.add_argument("--cort_yield_tens_func", type=str, default=None,
                        help="Function name for cortical tension yield. Defaults to --yield_tens_func.")
    parser.add_argument("--cort_poissons_ratio", type=float, default=None,
                        help="Poisson's ratio for cortical bone. Defaults to --poissons_ratio.")
    parser.add_argument("--pmma_E", type=float, default=2500,
                        help="Sets the Elastic Modulus for PMMA caps in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_v", type=float, default=0.3,
                        help="Sets the Poisson's ratio for the PMMA material(s) in the FE model. (default: %(default)s)")
    parser.add_argument("--pmma_thick", type=float, default=DEFAULT_PMMA_THICKNESS_MM,
                        help="Sets the minimum thickness for PMMA caps in the FE model. Default value (6 [mm]) based on observed measurement of Keaveny BCT FE modeling of the femur. (default: %(default)s [mm])")
    parser.add_argument("--pmma_intrusion", type=float, default=DEFAULT_PMMA_INTRUSION_MM,
                        help="Sets how far PMMA fixtures intrude toward the bone side before bone overlap is preserved during material combination. (default: %(default)s [mm])")
    parser.add_argument("--pmma_mat_id", type=int, default=5000,
                        help="Sets the material ID for the PMMA blocks. (default: %(default)s)")
    parser.add_argument("--fe_displacement", type=float, default=DEFAULT_FEMUR_FE_DISPLACEMENT,
                        help="Sets the applied displacement endpoint for the sideways-fall model. The default reports the force at 4%% displacement. (default: %(default)s)")
    parser.add_argument("--femur_shaft_length", type=float, default=DEFAULT_FEMUR_SHAFT_LENGTH_MM,
                        help="Retained proximal femur length [mm] for --femur_cut_mode fixed_length. (default: %(default)s [mm])")
    parser.add_argument("--femur_cut_mode", choices=["lesser_trochanter", "fixed_length"],
                        default=DEFAULT_FEMUR_CUT_MODE,
                        help="Set the distal femur crop. Default detects the lesser trochanter and fails if the field of view is incomplete. (default: %(default)s)")
    parser.add_argument("--femur_lesser_trochanter_distal_offset", type=float, default=DEFAULT_LESSER_TROCHANTER_DISTAL_OFFSET_MM,
                        help="Cut this many mm distal to the detected lesser trochanter in lesser_trochanter mode. (default: %(default)s [mm])")
    parser.add_argument("--femur_lesser_trochanter_distal_offset_percent", type=float, default=None,
                        help="Optional percentage of the detected greater-to-lesser trochanter z distance to cut distal to the lesser trochanter. Overrides --femur_lesser_trochanter_distal_offset when set. (default: %(default)s)")
    parser.add_argument("--femur_input_margin", type=float, default=DEFAULT_FEMUR_INPUT_MARGIN_MM,
                        help="Pad the input image/mask as needed so femur foreground has this margin before pre-rotation and ICP. (default: %(default)s [mm])")
    parser.add_argument("--compartment_mask", type=str, default=None,
                        help="Optional trabecular/cortical compartment mask aligned with the bone mask. Defaults: cortical=1, trabecular=2.")
    parser.add_argument("--cortical_label", type=int, default=DEFAULT_CORTICAL_LABEL,
                        help="Label value for cortical bone in --compartment_mask. (default: %(default)s)")
    parser.add_argument("--trabecular_label", type=int, default=DEFAULT_TRABECULAR_LABEL,
                        help="Label value for trabecular bone in --compartment_mask. (default: %(default)s)")
    parser.add_argument("--reference_path", type=str, required=False, default=None,
                        help="Path to the reference vtk file for ICP registration. (default: None)")
    parser.add_argument("--pmma_yield_compression", type=float, default=None,
                        help="Sets the yield strength in compression for PMMA material in the FE model. (default: %(default)s [MPa])")
    parser.add_argument("--pmma_yield_tension", type=float, default=None,
                        help="Sets the yield strength in tension for PMMA material in the FE model. (default: %(default)s [MPa])")
    # Parse and display
    args = parser.parse_args()
    
    # Set default reference paths
    if args.reference_path is None:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dat")
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
