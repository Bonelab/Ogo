# FEA code map

Start at `ogo/cli/GenerateFEM.py`. It parses `ogoFEA`, calls the anatomy builder,
audits the model, runs FAIM/N88 and writes modelling metadata.

| Task | File / function |
| --- | --- |
| CLI and site arguments | `GenerateFEM.build_parser`, `build_spine_command`, `build_femur_command` |
| Reporting endpoint and criteria | `GenerateFEM.solve_report_profile`, `critical_volume_percent`, `critical_strain` |
| Spine preparation and model | `spine._prepare_spine_inputs`, `_register_spine_inputs`, `process_vertebra` |
| Hip preparation and model | `femur._prepare_femur_inputs`, `_register_femur_inputs`, `sidewaysFallFe` |
| Hip multistart alignment | `hip_registration.estimate_femur_icp` |
| Physical scan-end/shaft measurements | `shaft_geometry.capture_distal_scan_face`, `measure_available_shaft`, `verify_model_shaft` |
| Shared transforms | `alignment.py`, `image_io.py` |
| PMMA/contact geometry | `boundary.generate_bone_cap_mask`, `generate_projected_material_disk_vtk` |
| Density laws and material tables | `material_laws.py`, `materials.build_bone_pmma_material_table` |
| N88 mesh and file writing | `model.create_microfe_model`, `write_model` |
| Solver discovery, Pistoia and results | `ogo/util/faim.py` |
| Generation measurements / QC rules | `validation.measure_model`, `evaluate` |
| Legacy measurement recovery | `legacy_qc.backfill_measurements` |
| Three-view previews | `qc_render.py`, using `ogo/cli/Visualize.py` |
| Gallery and review records | `gallery.py`, `review.py`, `gallery_review.js` |

Anatomy defaults are constants at the top of `spine.py` and `femur.py`.
Shared helpers own geometry/material logic; the CLI owns solving and reporting.
Changing a reporting criterion must not silently change model construction.
Hip supports only the GT-relative flat shaft crop; the initial fixed-length crop
is used for registration only. Spine settings are explicit keyword arguments;
misspelled setting names raise an error before reading inputs.

Spine cleans the body label, crops, estimates centroid-initialized ICP, resamples,
builds caps/materials and applies compression. Hip captures the original distal
scan face, rough-crops for registration, transforms the full scan, builds proximal
supports, checks available coverage, then crops and constrains the shaft.
Reference scaling changes the registration target, not native bone dimensions.
The modelling JSON records spine label smoothing and resolved contact settings.

For hip, the saved-model length check is
`min(z of Greater_Trochanter_PMMA_Nodes) - median(z of Distal_Femur_Nodes)`.
Coverage also accounts for the most proximal transformed scan-end point and
resampling clearance, so the final face must be complete, not an oblique remnant.

QC measurements are written during generation; validation reads these CSVs
without reopening models. Legacy recovery uses registered anatomy and saved
shaft geometry where available, recording what was recovered. Manual inclusion
is separate from the automatic construction checks.

Tests are in `tests/fea/` and `tests/util/`; run focused tests for the module
being edited. Imaging/mesh tests require vtkbone and SimpleITK. Document changes
to defaults and preserve their matching CLI/provenance tests.
