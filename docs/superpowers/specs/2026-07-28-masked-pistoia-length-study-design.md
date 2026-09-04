# Masked Pistoia and RETRO2 Femur Length Study Design

## Goal

Add a generic Pistoia region-of-interest mask workflow to Ogo so a solved FE
model can report both full-model Pistoia and masked-region Pistoia. The first
use case is the RETRO2 femur length study, where the mask is the femoral neck.
The same mechanism should also support spine models, for example vertebral-body
or trabecular-body-only Pistoia.

## Inputs

`ogo.cli.GenerateFEM` should accept a site-agnostic optional argument:

```text
--pistoia_mask PATH
```

The mask is interpreted in the input image space. For hip this can be a
femoral-neck mask. For spine this can be a vertebral-body, trabecular-body, or
other ROI mask. The mask is optional; if it is not provided, the current full
model Pistoia behavior is unchanged.

For the RETRO2 femur length study, inputs should come from the PI-preferred FEA
folder:

```text
/home/matthias.walle/data/scratch/RETRO2/fea/QCT/hip_QCT/
/home/matthias.walle/data/scratch/RETRO2/fea/segmentations/left_femur/
/home/matthias.walle/data/scratch/RETRO2/fea/segmentations/femoral_neck_original_dim/
```

Cases without a full hip QCT or without a reconstructed femoral-neck mask should
be reported in a manifest and skipped for masked femoral-neck Pistoia.

## Pipeline Behavior

The ROI mask must be carried through exactly the same spatial operations as the
FE model:

- input resampling
- rough crop or crop-to-mask steps
- padding
- ICP/reference transformation
- reference-grid resampling
- final length crop

ROI masks use nearest-neighbor interpolation. The transformed ROI mask should be
written beside the generated model for auditability. Its geometry must match the
final model grid so `ogo.util.faim.run_pistoia_with_image_mask()` can select
model elements directly.

The existing post-hoc source-space masks must not be applied directly to
post-ICP models. Prior checks showed that source-space masks select zero
elements in existing length-crop `.n88model` files.

## Pistoia Reporting

When `--run_pistoia` is enabled, Ogo should continue to run full-model Pistoia
with the existing critical strain and critical volume settings.

When `--pistoia_mask` is provided, Ogo should additionally run masked Pistoia
using the transformed model-space ROI mask. The lower-level mechanism should
reuse `ogo.util.faim.run_pistoia_with_image_mask()`, which temporarily excludes
elements outside the image mask and then restores the original model materials.

The output CSV should include both full and masked values. Generic column names
avoid femur-specific assumptions:

```text
pistoia_failure_load_N
pistoia_stiffness_N_per_mm
masked_pistoia_failure_load_N
masked_pistoia_stiffness_N_per_mm
masked_pistoia_critical_volume_pct
masked_pistoia_critical_ees
masked_pistoia_ees_at_crit_vol
masked_pistoia_mask_file
masked_pistoia_selected_elements
masked_pistoia_original_included_elements
```

If masked Pistoia fails and `--require_pistoia` is not set, the main solve
should still complete and report a warning, matching the current full-Pistoia
failure behavior.

## RETRO2 Length Study Output

Create a clean study folder:

```text
/home/matthias.walle/data/scratch/RETRO2/fea/length_study_femoral_neck_pistoia/
```

The folder should contain:

- `models/<subject>/<length>/...` for `.n88model` and sidecars
- `individual_results/<subject>/<length>/...` for solver and Pistoia outputs
- `transformed_masks/<subject>/<length>/...` for model-space Pistoia ROI masks
- `logs/` for generation and solve logs
- `tables/run_manifest.csv`
- `tables/results_summary.csv`

`results_summary.csv` should retain the length-study fields from the existing
study and add the masked Pistoia fields above.

## Validation

Implementation should include unit or integration tests that verify:

- `--pistoia_mask` appears in the hip and spine CLI.
- The femur pipeline resamples/transforms the ROI mask to the final model grid.
- The transformed ROI mask selects nonzero model elements in a generated model.
- Full Pistoia and masked Pistoia can both be reported in the same results CSV.
- Missing ROI masks produce manifest warnings rather than silent omission.

Before running the full RETRO2 study, run a small smoke set with available
femoral-neck masks across several crop lengths and verify:

- transformed ROI masks overlay the final model anatomy;
- selected-element counts are plausible and nonzero;
- full-model and masked Pistoia outputs are both present.
