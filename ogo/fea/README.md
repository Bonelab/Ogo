# Ogo FEA Workflow

This directory contains the maintained model-building code used by
`ogoFEA`. The user-facing entry point is the CLI in
`ogo/cli/GenerateFEM.py`; the lower-level spine and femur modules in this
directory should stay implementation details.

The workflow builders are organized by anatomy:

- `spine.py` owns spine-compression defaults, preprocessing, model building,
  and the lower-level spine command parser.
- `femur.py` owns hip sideways-fall defaults, preprocessing, model building,
  and the lower-level femur command parser.
- Shared helper modules such as `alignment.py`, `boundary.py`, `materials.py`,
  and `model.py` contain reusable mechanics used by both anatomy workflows.

## Inputs

All FEA workflows expect a density-calibrated CT image and an aligned
segmentation image.

The calibrated image should be a scalar image in K2HPO4-equivalent density
units, with the same image space as the segmentation. During model generation,
the image is cropped, aligned to the model reference frame, resampled to
isotropic spacing, thresholded at -31 mg/cc, and converted to material IDs.

The segmentation depends on the model type:

| Model | Required mask | Label convention |
| --- | --- | --- |
| Spine compression | Labelled vertebra mask | Provide the vertebral body label and posterior process label for each target as `LEVEL:BODY_LABEL:PROCESS_LABEL`. |
| Hip sideways fall | Whole-femur mask | Nonzero voxels are the femur by default. The wrapper chooses left/right from `--side`. |
| Hip with trab/cort compartments | Whole-femur mask plus `--compartment_mask` | Compartment mask defaults are cortical `1`, trabecular `2`; override with `--cortical_label` and `--trabecular_label` if needed. |

The spine mask must contain both the body and posterior process labels. The
femur workflow currently requires a side-specific run because alignment uses a
side-specific femur reference frame.

## Implementation Walkthrough

The high-level command is `ogoFEA`, implemented in `ogo/cli/GenerateFEM.py`.
That file is intentionally the best starting point for a new developer because
it owns the user-facing arguments, the lower-level anatomy command construction,
the solve/reporting path, and the `_modeling.json` provenance record.

Use this map when tracing a run:

| Question | Main code location |
| --- | --- |
| Which command-line options are exposed? | `ogo/cli/GenerateFEM.py::build_parser` |
| How does `ogoFEA spine` become a lower-level model-builder call? | `ogo/cli/GenerateFEM.py::build_spine_command` |
| How does `ogoFEA hip` become a lower-level model-builder call? | `ogo/cli/GenerateFEM.py::build_femur_command` |
| Where are solve endpoints, Pistoia settings, and result columns selected? | `ogo/cli/GenerateFEM.py::solve_report_profile`, `critical_volume_percent`, `solve_model` |
| Where is model provenance written? | `ogo/cli/GenerateFEM.py::write_modeling_metadata` |
| Where is the generated `.n88model` audited? | `ogo/cli/GenerateFEM.py::audit_generated_model`, `ogo/cli/CheckFEModelBC.py` |

The anatomy-specific implementation is split below the CLI:

| Anatomy | Main file | What belongs there |
| --- | --- | --- |
| Spine compression | `ogo/fea/spine.py` | Vertebra labels, body/process QC, spine ICP, 1 mm resampling, stable PMMA cap generation, axial compression defaults. |
| Hip sideways fall | `ogo/fea/femur.py` | Side handling, femur reference alignment, rough pre-ICP crop, post-ICP GT-length crop, femoral-head/GT/distal supports, sideways-fall defaults. |
| Shared material tables | `ogo/fea/materials.py` | Bone material ID ranges, PMMA material, shared femur/spine material-table construction. |
| Material laws | `ogo/fea/material_laws.py` | Density-to-modulus and yield-strength functions such as `default_E` and `kopperdahl_trab_E`. |
| Shared geometry/BC helpers | `ogo/fea/boundary.py` | Resampling helpers, projected material disks, PMMA caps, contact footprint calculations. |
| N88 model writing | `ogo/fea/model.py` | `create_microfe_model` and `write_model`. |
| FAIM/N88 adapter | `ogo/util/faim.py` | Solver command resolution, running `faim`/N88 tools, Pistoia parsing, masked Pistoia, compact `_results.csv`. |

For spine, the maintained path is:

1. `GenerateFEM.py::build_spine_command` forwards the calibrated image,
   segmentation, target `LEVEL:BODY_LABEL:PROCESS_LABEL`, preset, optional
   Pistoia mask, and any lower-level overrides.
2. `spine.py::main` thresholds the requested body/process labels, crops around
   the vertebra, checks posterior-process orientation, and runs scaled ICP to
   the bundled L4 body reference from `default_spine_reference_path`.
3. `spine.py` applies the transform and resamples the density image with cubic
   interpolation and labels/masks with nearest-neighbor interpolation.
4. `boundary.py::generate_bone_cap_mask` and related helpers generate superior
   and inferior PMMA caps. The stable-contact settings are defined by
   `DEFAULT_SPINE_STABLE_CONTACT_*` constants in `spine.py`.
5. `materials.py::build_spine_material_table` assigns bone and PMMA material
   definitions. The default benchmark linear preset uses
   `material_laws.kopperdahl_trab_E`.
6. `model.py::create_microfe_model` writes the N88 model structure, and the
   high-level wrapper solves at the selected endpoint, defaulting to Crawford
   `0.68%` strain for the maintained linear spine workflow.

For hip, the maintained short-femur path is:

1. `GenerateFEM.py::build_femur_command` forwards side, optional compartment
   mask, optional Pistoia mask, and the femur crop settings.
2. `femur.py::main` selects the side-specific femur mask and aligns it to the
   bundled left or mirrored-right reference. The fixed rough crop for ICP
   stability is controlled by `DEFAULT_FEMUR_ROUGH_PRE_ICP_LENGTH_MM`.
3. In `greater_trochanter_length` mode, the rough crop is used only for the ICP
   estimate. The final model is then regenerated from the transformed full scan.
4. `femur.py::crop_vtk_images_to_greater_trochanter_length` applies the final
   post-ICP flat shaft crop. The requested distal length is measured after the
   generated GT support, not from the femoral head and not from the raw GT
   landmark.
5. `femur.py::proximal_sideways_fall_fixture_plane` and
   `boundary.py::generate_projected_material_disk_vtk` create the femoral-head
   and greater-trochanter PMMA supports. The distal crop face becomes the distal
   shaft support.
6. `materials.py::build_femur_material_table` assigns bone and PMMA material
   definitions. The simple whole-femur path uses one trabecular-style bone
   region unless a compartment mask is supplied.
7. The high-level wrapper solves and reports the maintained hip endpoint,
   defaulting to `4%` displacement.

The most important short-femur convention is the final shaft-length definition.
In `greater_trochanter_length` mode the requested value is the distance between
the final GT support and the final distal shaft boundary condition. The audit in
`ogo/cli/CheckFEModelBC.py::audit_femur_sideways` reports the robust model-space
measurement as:

```text
p5(z of Greater_Trochanter_PMMA_Nodes) - median(z of Distal_Femur_Nodes)
```

This is deliberately a boundary-condition measurement. It verifies the
mechanical model that FAIM/N88 will solve rather than an intermediate crop
coordinate.

## Basic Commands

### Recommended Spine Command

For a standard solved L1 compression model, run the high-level wrapper. This
generates the `.n88model`, audits the boundary conditions, solves it with
FAIM/N88, runs optional Pistoia postprocessing, and writes a compact
`_results.csv` next to the model:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --output_path derivatives/fea \
  --threads 4 \
  --faim_bin_dir /path/to/n88/bin \
  --run_pistoia \
  --critical_volume 12 \
  --critical_strain 0.007
```

The `--vertebra` value is `LEVEL:BODY_LABEL:PROCESS_LABEL`. For example,
`L1:2:3` means that label `2` is the L1 vertebral body and label `3` is the L1
posterior process. Repeat `--vertebra` to process more than one level from the
same image.

To additionally report Pistoia for an ROI such as the vertebral body, pass a
mask in the same image space as the calibrated CT and segmentation:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --pistoia_mask sub-001_desc-L1Body_mask.nii.gz \
  --output_path derivatives/fea \
  --threads 4 \
  --faim_bin_dir /path/to/n88/bin \
  --run_pistoia \
  --critical_volume 12 \
  --critical_strain 0.007
```

Ogo transforms this ROI mask through the same preprocessing, alignment,
resampling, and padding steps as the model, then writes the model-space mask as
`<model-stem>_pistoia_mask.nii.gz`. The `_results.csv` includes both full-model
Pistoia fields and `masked_pistoia_*` fields when a mask is supplied.

Generate and solve one L1 model:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --output_path derivatives/fea
```

Generate L1 and L2 in one command:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --vertebra L2:4:5 \
  --output_path derivatives/fea
```

Generate only the model files, without running FAIM:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --output_path derivatives/fea \
  --no-solve
```

Generate the left femur:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-leftFemur_mask.nii.gz \
  --side left \
  --output_path derivatives/fea
```

Generate both femurs from a mask containing both sides:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-femur_mask.nii.gz \
  --side both \
  --output_path derivatives/fea
```

Generate a femur model with an explicit trabecular/cortical compartment mask:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-leftFemur_mask.nii.gz \
  --side left \
  --compartment_mask sub-001_desc-leftFemur_compartments.nii.gz \
  --output_path derivatives/fea
```

Use `--dry-run` to print the lower-level command without generating files.

## Solver Settings

By default, `ogoFEA` generates each `.n88model`, audits the boundary
conditions, writes a `_modeling.json` record, then runs the FAIM/N88 solve and
postprocessing.

Control the FAIM thread count with `--threads`:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --threads 8
```

This sets the FAIM solver command to use the multithreaded engine with
`--threads=8`. It is a per-solve thread count, not a request to run multiple
models in parallel.

If Ogo and FAIM/N88 are installed in separate conda environments, activate the
Ogo environment first and pass the FAIM environment with `--faim_env`. Ogo
generates the `.n88model` in the current Python process, then launches each
FAIM/N88 subprocess through `conda run -n ENV_NAME`.

For example, if Ogo is installed in `ogo-dev` and FAIM/N88 is installed in
`ogoloco-n88`:

```bash
conda activate ogo-dev

ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --faim_env ogoloco-n88 \
  --threads 8
```

The same pattern works for femur:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-leftFemur_mask.nii.gz \
  --side left \
  --faim_env ogoloco-n88 \
  --threads 8
```

Use `--conda_executable` if `conda` is not on `PATH` or if a specific conda
install should be used. For example:

```bash
ogoFEA spine \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-spine_labels.nii.gz \
  --vertebra L1:2:3 \
  --faim_env ogoloco-n88 \
  --conda_executable /path/to/conda
```

Other available solver-location options are `--faim_bin_dir`,
`--faim_install_root`, `--faim_license_dir`, and explicit command overrides
such as `--faim_command`, `--n88postfaim_command`, and `--n88tabulate_command`.
Ogo does not search your filesystem for an arbitrary N88 installation. If the
N88 commands are not already on `PATH`, use `--faim_bin_dir`, `--faim_env`, or
the explicit command overrides.

FAIM/N88 does not interpret prescribed displacement as percent strain in these
models. The `.n88model` stores the applied displacement in millimeters. For
spine, Ogo converts Crawford `0.68%` strain to millimeters from the generated
vertebral height before solving, updates the model displacement, and reports
that solved load directly. Femur reporting still converts the `4%` endpoint
from generated model geometry for `_results.csv`.

## Outputs

The default output directory is the calibrated image directory. Use
`--output_path` for a dedicated derivatives folder.

Expected default outputs are:

| Output | Meaning |
| --- | --- |
| `.n88model` | The generated finite-element model. |
| `_modeling.json` | Traceable record of inputs, alignment, image processing, material laws, boundary conditions, audit summary, and reporting settings. |
| `_results.csv` | One-row summary written after solve with stiffness, reaction force, and profile-specific force endpoints. |
| Solve/postprocessing files | Written when `--no-solve` is not used. These come from FAIM/N88 postprocessing. |

Use `--debug` to write visual sidecars such as the boundary-condition audit
PNG. These images are intended for manual model inspection and are not part of
the minimal output contract.

## Results CSV

After a successful solve, `ogoFEA` writes one compact `_results.csv` next to the
model. This is the file intended for downstream analysis tables.

Maintained spine and femur models write the same columns:

| Field | Meaning |
| --- | --- |
| `model_file` | Source `.n88model`. |
| `analysis_file` | FAIM/N88 postprocessing text file used to tabulate reaction force. |
| `analysis_var` | N88 reaction-force variable tabulated from `analysis_file`. |
| `applied_displacement` | Prescribed displacement in the solved `.n88model`. |
| `reaction_force_N` | Reaction force from the solved model at the prescribed displacement. |
| `stiffness_N_per_mm` | Reaction force divided by applied displacement when computable. |
| `characteristic_length_mm` | Spine: distance between `body_top` and `body_bottom` centroids along the loading axis. Femur: full generated-model span along the loading axis. |

## Reaction Force and Displacement Endpoint

The reported `reaction_force_N` is the net support reaction from the solved
model at the prescribed displacement endpoint. For spine compression this is
tabulated from the axial `fz_ns1` reaction force. For hip sideways fall this is
tabulated from the lateral `fy_ns1` reaction force. The sign depends on the
model coordinate convention, so the compact CSV reports the force magnitude
used for stiffness and downstream comparisons.

The reaction force and stiffness in `_results.csv` should be interpreted at
that endpoint. They are the relevant summary values for the chosen model
protocol: the force is the load carried by the model at the target normalized
displacement, and stiffness is the secant stiffness from zero displacement to
that same endpoint. Changing the endpoint can change both values, especially
for nonlinear material settings, so comparisons should use the same model type,
material law, boundary conditions, and displacement endpoint.

Ogo chooses percent displacement endpoints so models of different physical
size are compared at a similar normalized deformation:

| Model | Default endpoint | Length used for conversion | Reason |
| --- | --- | --- | --- |
| Spine compression | `0.68%` strain | Distance between `body_top` and `body_bottom` centroids along the loading axis | Matches the vertebral strength endpoint reported by Crawford et al. (2003), and keeps the solved displacement scaled to each vertebral body height. |
| Hip sideways fall | `4%` displacement | Full generated-model span along the loading axis | Keeps the sideways-fall load endpoint scaled to each generated femur model instead of using one fixed millimeter displacement for all femora. |

FAIM/N88 stores prescribed displacement in millimeters, not percent strain.
Before solving, `ogoFEA` converts the percent endpoint to millimeters from the
generated model geometry and updates the `.n88model`. The exact percent
endpoint, converted millimeter displacement, and characteristic length are
recorded in `_modeling.json`; `_results.csv` stays compact for downstream
statistics. `stiffness_N_per_mm` is computed as
`reaction_force_N / applied_displacement`.

## Reference Models

The reference models define the common coordinate frame used before meshing and
boundary-condition placement. They are alignment targets, not subject-specific
geometry copied into the final model.

For spine compression, the bundled vertebral-body reference is an L4 body. The
same reference is used for all requested vertebrae. Before ICP, Ogo estimates
principal-axis lengths for the segmented body and for the L4 reference, scales
the reference to match the target vertebra within the configured scale limits,
and then runs rigid ICP. This lets one L4 reference provide a stable model frame
for other vertebral levels while still adapting to their size.

For hip sideways fall, the canonical bundled reference is a left femur. Left
femur models ICP-align directly to that reference. Right femur models use
`RT_FEMUR_SIDEWAYS_FALL_REF.vtk` when it is present; otherwise Ogo mirrors the
left femur reference across the x direction in memory and uses that mirrored
surface as the right-side ICP target. The mirrored reference keeps left and
right outputs in matching sideways-fall coordinate frames without requiring a
separate right reference file.

## Spine Model

The spine workflow builds one compression model per `--vertebra` target:

1. Threshold the labelled mask into vertebral body and posterior process.
2. Crop the image and masks to the vertebra.
3. Align the body to the scaled L4 vertebral-body reference using ICP.
4. Apply transform and isotropic resampling through the shared VTK reslice
   helper. The default output spacing is `1.0 x 1.0 x 1.0 mm`.
5. Smooth body/process masks with one binary close/open pass only when at
   least one input spacing dimension is coarser than 2 mm.
6. Generate superior and inferior PMMA caps fitted to the body surface.
7. Convert density to material IDs and construct the N88 model.
8. Apply axial compression boundary conditions.

Default spine boundary conditions:

| Region | Node set | Constraint |
| --- | --- | --- |
| Superior PMMA cap | `body_top` | Prescribed displacement in the compression direction. |
| Inferior PMMA cap | `body_bottom` | Fixed in all displacement directions. |

The default high-level spine preset is `benchmark-linear`. It applies the
spineFE benchmark linear material settings for model construction, then
`ogoFEA` updates the generated model to solve at Crawford `0.68%` strain.
The `0.68%` endpoint follows Crawford RP, Cann CE, Keaveny TM. Finite element
models predict in vitro vertebral body compressive strength better than
quantitative computed tomography. Bone. 2003;33(4):744-50.
doi:10.1016/s8756-3282(03)00210-2.

Available spine presets:

| Preset | Use |
| --- | --- |
| `benchmark-linear` | Default linear spineFE benchmark-style run. |
| `benchmark-nonlinear` | Nonlinear material/yield settings for benchmark-style runs. |
| `none` | Pass only explicit lower-level options. |

## Femur Model

The femur workflow builds one sideways-fall model per side:

1. Read the calibrated image and whole-femur mask.
2. Pre-rotate the side to a stable starting orientation.
3. ICP-align to the left femur reference or the mirrored right reference.
4. Apply transform and isotropic resampling through the shared VTK reslice
   helper. The default output spacing is `1.0 x 1.0 x 1.0 mm`.
5. Smooth the transformed femur mask with one binary close/open pass only when
   at least one input spacing dimension is coarser than 2 mm. If a compartment
   mask is supplied, the derived cortical binary mask follows the same rule.
6. Standardize the distal shaft with a flat post-ICP model-grid cut. In
   `greater_trochanter_length` mode, the fixed rough crop is used only to make
   ICP stable; the final model is then generated from the transformed full scan.
   The requested shaft length is measured after the generated
   greater-trochanter support, not from the femoral-head support and not from
   the raw GT landmark alone.
7. Generate PMMA fixtures for the femoral head and greater trochanter. The
   femoral-head fixture is placed on the high-y side and the greater-trochanter
   fixture is placed on the low-y side. Both fixtures are anchored to the
   proximal model footprint so they do not overwrite bone voxels.
8. Convert density to material IDs and construct the N88 model.
9. Apply sideways-fall boundary conditions.

Default femur boundary conditions:

| Region | Node set | Constraint |
| --- | --- | --- |
| Femoral head PMMA | `Femoral_Head_PMMA_Nodes` | Prescribed displacement toward the greater trochanter. |
| Greater trochanter PMMA | `Greater_Trochanter_PMMA_Nodes` | Constrained in the loading direction. |
| Distal shaft cut face | `Distal_Femur_Nodes` | Constrained to remove rigid-body motion. |

Femur reporting defaults to stiffness and force at `4%` displacement. When
`--run_pistoia` is supplied, the high-level wrapper also reports full-model
Pistoia failure load and, if `--pistoia_mask` is supplied, masked-region
Pistoia failure load.

For the maintained short-femur cohort model, use:

```bash
ogoFEA hip image.nii.gz femur_mask.nii.gz \
  --side left \
  --femur_cut_mode greater_trochanter_length \
  --femur_greater_trochanter_distal_length 20 \
  --femur_shaft_length 120 \
  --run_pistoia \
  --critical_volume 10 \
  --critical_strain 0.009
```

The generated model can be audited with `ogoFEA check-bc`. For hip sideways
fall, the post-GT-support shaft length is measured from final boundary-condition
node sets as:

```text
p5(z of Greater_Trochanter_PMMA_Nodes) - median(z of Distal_Femur_Nodes)
```

The low GT percentile and distal-shaft median make this measurement robust to a
small number of edge voxels or support-surface bleed-through nodes.

## Materials

Spine and femur use the same shared material-table builder in
`ogo.fea.materials`.

Material ID convention:

| Material | IDs |
| --- | --- |
| Background | `0` |
| Trabecular bone | `1..128` |
| Cortical bone | `129..256` when a cortical region is present |
| PMMA fixtures/caps | Default `5000` |

The default material law is `default_E`, defined in `ogo.fea.material_laws`:

```text
E = 10500 * rho^2.29 MPa
```

where `rho` is density in g/cc. Poisson's ratio defaults to `0.3` for bone.
PMMA defaults are `E = 2500 MPa`, Poisson's ratio `0.3`.

The default spine preset overrides the elastic law to `kopperdahl_trab_E` for
both trabecular and cortical regions. The simple femur path uses `default_E`
unless a material-law override is supplied. If a femur compartment mask is
provided, cortical and trabecular regions use the same law unless cortical
overrides are supplied.

### Kopperdahl Spine Presets

The spine benchmark presets use the Kopperdahl trabecular density-modulus law:

```text
E = 2980 * rho^1.05 MPa
```

where `rho` is density in g/cc. In code this is
`ogo.fea.material_laws.kopperdahl_trab_E`.

The `benchmark-linear` preset uses this elastic law only. It does not assign a
bone yield criterion, so the solved model is linear elastic. The maintained
spine workflow runs Pistoia when `--run_pistoia`, `--require_pistoia`, or
`--pistoia_mask` is supplied.

The `benchmark-nonlinear` preset additionally assigns the Kopperdahl
compressive yield law:

```text
sigma_y,c = 37.4 * rho^1.39 MPa
```

where `rho` is density in g/cc. In code this is
`ogo.fea.material_laws.kopperdahl_trab_yc`. The current nonlinear preset uses
the same function for compression and tension unless the user supplies separate
yield functions. When yield functions are present, Ogo creates
`vtkboneMohrCoulombIsotropicMaterial` entries with the resolved tensile and
compressive yield strengths.

Material-law overrides can be passed through the high-level command because
unknown options are forwarded to the lower-level generator:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-leftFemur_mask.nii.gz \
  --side left \
  --elastic_E_func kopperdahl_trab_E \
  --cort_elastic_E_func kopperdahl_trab_E
```

Supported named laws are defined in `ogo.fea.material_laws`, including
`default_E`, `kopperdahl_trab_E`, `kopperdahl_trab_yc`, `morgan_trab_E`,
`crawford_voxel_E`, and Bayraktar-style constants.

## Quality Control

Every generated model is audited by `ogo/cli/CheckFEModelBC.py` unless
`--skip_bc_audit` is used. The audit summary is embedded in `_modeling.json`.

For manual inspection, run with `--debug`:

```bash
ogoFEA hip \
  sub-001_desc-vqct_ct.nii.gz \
  sub-001_desc-leftFemur_mask.nii.gz \
  --side left \
  --debug \
  --no-solve
```

This writes a boundary-condition audit PNG next to the model. For spine, debug
mode also enables quick-look output from the lower-level generator.

## Common Adjustments

Change output resolution:

```bash
ogoFEA spine image.nii.gz labels.nii.gz \
  --vertebra L1:2:3 \
  --iso_resolution 0.8
```

Change femur shaft retention after the greater-trochanter support:

```bash
ogoFEA hip image.nii.gz femur_mask.nii.gz \
  --side left \
  --femur_cut_mode greater_trochanter_length \
  --femur_greater_trochanter_distal_length 45
```

Change PMMA dimensions:

```bash
ogoFEA hip image.nii.gz femur_mask.nii.gz \
  --side left \
  --pmma_thick 6 \
  --pmma_intrusion 6
```

For hip sideways fall, both PMMA fixtures use bbox-relative contact planes by
default. The maintained footprint is rectangular, centered from stored bbox
fractions, and scaled independently to `1.1 x 1.1` of the model bbox in the two
authored in-plane axes. `--pmma_thick` controls total fixture thickness;
`--pmma_intrusion` controls how far anatomy can occupy that fixed thickness
before unsupported contact columns are removed. The PMMA labels do not replace
bone voxels.

For spine compression, `--pmma_thick` controls the total generated disk
thickness and `--pmma_intrusion` controls how far anatomy can occupy that fixed
thickness. The maintained spine defaults are 10 mm thickness and 6 mm intrusion.
Intrusion does not overwrite vertebral-body voxels; it limits how far from the
superior/inferior body surface the cap is allowed to search for supporting
anatomy before generating the flat PMMA disk.

After vtkbone identifies visible cap nodes, spine compression filters each cap
node set to the dominant coordinate plane along the load axis. This removes
small rim/contact-side nodes from the jelly-bean disk surface while keeping the
full disk material geometry in the solved model.

```bash
ogoFEA spine image.nii.gz labels.nii.gz \
  --vertebra L1:20:48 \
  --pmma_thick 10 \
  --pmma_intrusion 6
```

Run model generation only, inspect the `.n88model`, then solve later:

```bash
ogoFEA hip image.nii.gz femur_mask.nii.gz \
  --side left \
  --no-solve \
  --debug
```
