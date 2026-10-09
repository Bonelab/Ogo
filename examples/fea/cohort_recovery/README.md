# Spine Cohort Recovery

This study-specific example resegments selected cases on full native HU CTs,
isolates L1 (vertebral-level label 20), and maps body 2 / process 1 into the
calibrated QCT's physical grid. Only then are inputs cropped with 15 mm padding.
It preserves original study outputs and builds a separate after-only review ZIP.
Successful solving is not evidence of anatomical validity: review the new models.

Use a GPU compute node, never an HPC login node. Install Ogo and spine-segment
separately and use the validated segmentation model bundle for your study.

```bash
python examples/fea/cohort_recovery/recover_spines.py \
  --root /scratch/study/spine_recovery \
  --review-csv /scratch/study/study_inclusion.csv \
  --selection manual-exclusions \
  --library /scratch/study/ct-library \
  --qct-dir /scratch/study/spine_QCT \
  --segment-cli /path/to/spine-segment \
  --model-bundle /path/to/model-bundle-pytorch \
  --faim-bin-dir /path/to/n88/bin \
  --gpu-workers 2 --fea-workers 4
```

Use `--selection incomplete` for `incomplete_fea` cases not already manually
excluded. Use `--preflight` to check selected inputs without starting inference.
Native layout: `sub-ID/ses-1/ct/sub-ID_ses-1_ct.nii.gz`; QCT: `ID_QCT.nii.gz`;
review model ID: `ID_QCT_L1`. This example is not a generic BIDS discovery tool.

Defaults reproduce the study recovery: four CPU threads per solve, full-bone
critical volume 12%, body critical volume 35%, critical strain 0.007. Registration
and supports use the checked-out Ogo defaults. The launch record saves the Ogo
commit and segmentation bundle manifest hash. Preserve the bundle and its source
release separately. Do not change code or bundle mid-batch.

Per-case status and logs support resuming the same command. Use a new output root
when changing the selected cases or processing settings. Existing completed stages
are retained. Failures stay visible in `after_review.zip`; old manual decisions
are provenance only, not new inclusion decisions. `models/` contains n88models,
individual result and measurement CSVs, and before/after-solve WebP previews.

Participant lists, trained weights, scans, and study results are not included in
this repository. Local delivery/download utilities are operational helpers rather
than part of the scientific model-generation pipeline.
