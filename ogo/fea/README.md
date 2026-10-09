# Hip and spine FEA

`ogoFEA` generates, checks and solves N88 models from calibrated QCT and
segmentations. Start with the [installation and L1 example](../../docs/fea/vertebra_fea_quickstart.md).
For implementation details, see the [code map](../../docs/fea/implementation.md).

## Run a model

Inputs are NIfTI images in the same physical space. QCT density must be in
mg/cm3 K2HPO4-equivalent units, not HU. Spine labels identify the body and
posterior process separately; a hip mask contains the selected whole femur.
Label numbers are dataset-specific.

```bash
# L1: body label 20, process label 48; also evaluate body-only Pistoia.
ogoFEA spine qct.nii.gz labels.nii.gz --vertebra L1:20:48 \
  --pistoia_mask_label 20 --critical_volume 12 --masked_critical_volume 35 \
  --critical_strain 0.007 --threads 4 --output_path derivatives/fea

# Left femur with a separate femoral-neck mask.
ogoFEA hip qct.nii.gz left_femur.nii.gz --side left \
  --pistoia_mask femoral_neck.nii.gz --critical_volume 11.2 \
  --masked_critical_volume 35 --critical_strain 0.009 \
  --threads 4 --output_path derivatives/fea
```

A regional mask enables both full-bone and masked Pistoia. Without a mask,
add `--run_pistoia` for full-bone results. `--require_pistoia` makes a failed
Pistoia calculation fatal. `--pistoia_mask_label` without a separate mask
selects labels from the input segmentation; repeat it for multiple labels.
Without label selection, all nonzero mask voxels are included.

Useful options:

- `--no-solve`: generate the model and QC without FAIM.
- `--dry-run`: print generator/solver commands without writing models.
- `--vertebra`: repeat for additional levels; format `LEVEL:BODY:PROCESS`.
- `--side`: `left`, `right` or `both` (default both); choose explicitly.
- `--threads`: per-job generation and solver threads, not parallel cases.

Use `ogoFEA spine --help` or `ogoFEA hip --help` for wrapper options.
Additional options are forwarded to the anatomy builder; its help is available
with `python -m ogo.fea.spine --help` or `python -m ogo.fea.femur --help`.

## Solver setup

FAIM/N88 tools must be on `PATH`, or supplied through `--faim_bin_dir DIR`,
`--faim_install_root DIR` or `--faim_env ENV`. Individual command overrides
are also available. Ogo does not search the filesystem for installations.
Use `--faim_license_dir DIR` when the license is not already configured.
On ARC, `source /work/boyd_lab/conda_init/ogo_fea.sh` configures the shared
environment; run generation and solving on a compute node.

## Model conventions

Both workflows use 1-mm voxels, cubic density interpolation and nearest-neighbor
labels. Reference scaling estimates registration; the solved subject geometry
is transformed rigidly, not scaled.

| Setting | Hip sideways fall | Spine compression |
| --- | --- | --- |
| Registration | NumPy ICP, 40,000 points, four axial starts | Centroid-initialized NumPy ICP, 8,000 points, 50 iterations |
| Target | Side-specific femur surface | L4 body surface |
| PMMA thickness / intrusion | 10 / 6 mm | 10 / 6 mm |
| Bone / PMMA Poisson ratio | 0.3 / 0.3 | 0.3 / 0.3 |
| PMMA modulus | 2500 MPa | 2500 MPa |
| Reporting endpoint | 4% of model span along loading axis | 0.68% of top-to-bottom BC-centroid distance |
| Full-bone Pistoia default | 11.2% volume, EES 0.009 | 2% volume, EES 0.007 |

The spine command above explicitly uses 12%, rather than the generic 2%
fallback. Regional critical volume inherits the full-bone value unless set
with `--masked_critical_volume`; it is not automatically 35%.
Critical volume is a percentage; critical strain is a dimensionless fraction.
These are analysis choices, not universally validated failure thresholds.

Hip registration uses a 120-mm proximal rough crop, then returns to the
transformed full scan. Proximal supports are built before the flat distal crop.
The default shaft is 10 mm distal to the generated GT support edge, after
clearing the transformed incomplete scan end. Final length is verified from
GT and distal BC nodes. Insufficient coverage is recorded as `too_short` and
is not solved. `--femur_greater_trochanter_distal_length` changes this length;
the 11.2% criterion is specific to the 10-mm protocol, not other lengths.

Spine keeps the largest face-connected body component without changing process
labels. Registration rejects implausible axial/process rotations rather than
returning a least-bad fit. Caps use axial body contact, small-tip trimming,
gap closing and largest-component retention without overwriting bone. Default
trim fraction is 0.10, minimum shift 3 mm. See `boundary.py` and the
`DEFAULT_SPINE_STABLE_CONTACT_*` constants for the exact construction.
Optional `--registration_backend vtk` and `--spine_icp_target vertebra` are
available; whole-vertebra ICP requires a matching aligned reference surface.

For linear materials, the hip converts QCT to ash density
`rho_ash = 1.06 * rho_QCT + 0.0389` (g/cm3), then uses
`E = 10500 * rho_ash^2.29` MPa. The default spine `benchmark-linear` preset
uses `E = 2980 * rho_QCT^1.05` MPa without ash conversion. Named laws are in
`material_laws.py`; `--elastic_E_func` overrides the law. Spine also supports
`benchmark-nonlinear` (yield materials, 4% endpoint) and `none` presets;
do not mix these with linear reference results.

## Outputs and QC

| Suffix | Contents |
| --- | --- |
| `.n88model` | Generated model, updated with solution fields after solving |
| `_results.csv` | Reaction force, stiffness and requested Pistoia results |
| `_modeling.json` | Inputs, processing, materials, BCs and solver settings |
| `_qc_metrics.csv` | Generation measurements, extended with solved SED checks |
| `_qc_3d.webp`, `_sed_3d.webp` | Before/after three-view previews |
| `_shaft_geometry.json` | Hip available, requested and verified retained length |
| `_pistoia_mask.nii.gz` | Regional mask transformed into model coordinates |

CSV forces are in N and stiffness in N/mm; Pistoia also has kN fields.
Reaction-force sign follows the model convention; use magnitudes for comparisons.
Linear stiffness is reaction-force magnitude divided by applied displacement;
Pistoia stiffness is reported separately. Missing results are blank, not zero.

Previews use white backgrounds, gold supports and red constrained surfaces.
Solved bone uses a shared linear Jet SED scale. Display smoothing does not
change FE results. A graphics failure warns but does not invalidate the model.
Compact WebP previews are for QC, not quantitative analysis. Legacy PNGs are
supported. Too-short hips receive an available-anatomy preview, not a solved model.

```bash
ogoValidateFEA derivatives/fea --output derivatives/qc --gallery-zip --workers 4
```

This writes `qc_summary.csv`, `study_inclusion.csv` and portable `gallery.zip`.
Extract the ZIP and open `gallery.html`; models are not bundled. Use `--gallery`
instead for an unpacked gallery and `--reviews study_inclusion.csv` to restore
manual decisions. Filters, sorting, page size and individual before/after views
support large-cohort review. Include/exclude records a reason and reviewer.
Export the reviewed CSV: browser edits do not update the project CSV or ZIP.

QC fails construction errors, flags uncertain anatomy or missing evidence for
review, and does not apply strength thresholds. A pass is not proof of correct
segmentation. Hip tests include shaft length, complete distal section, intended
support-patch coverage and X/Z fixation; spine tests include cap contact and
body/process orientation. Substantial body cleanup is flagged for review.

Old models without measurements are backfilled once and cached. Use `--site`
when anatomy cannot be inferred and `--evidence evidence.csv` for registered
body/process masks or shaft sidecars stored elsewhere. Evidence columns are
`site,model_id,body_mask,process_mask,shaft_geometry`; paths are relative to the
CSV. Native-space masks cannot substitute for model-space evidence. Missing
scan coverage stays unknown; backfill neither rewrites models nor renders images.
Failures before meshing require job logs for participant accounting.

For full-size resegmentation of reviewed failures, see the
[cohort recovery example](../../examples/fea/cohort_recovery/README.md).
