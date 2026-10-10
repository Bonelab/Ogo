# L1 FEA quickstart

## Install

Install Git and conda first, then:

```bash
git clone https://github.com/Bonelab/Ogo.git
cd Ogo
git checkout ogo-fea
conda create -n ogo -c numerics88 -c conda-forge \
  git pandas nibabel pyyaml vtkbone netcdf4 pbr "setuptools<69" wheel python=3
conda activate ogo
PBR_VERSION=0.0.1 python -m pip install -e .
ogoFEA --help
```

`PBR_VERSION` supplies installation metadata only. The Numerics88 channel
provides vtkbone; solving additionally requires FAIM/N88 tools and a license.
Use tools already on `PATH`, or add `--faim_bin_dir /path/to/n88/bin` to commands.
A separate solver conda environment can be selected with `--faim_env NAME`.

## Generate and inspect

Use density-calibrated QCT in mg/cm3 K2HPO4-equivalent units and a segmentation
in the same physical space. This example uses body label 20 and process label 48;
replace these with your dataset's labels. NIfTI geometry must preserve anatomical
orientation, with physical z along the superior-inferior axis.

```bash
ogoFEA spine qct.nii.gz labels.nii.gz --vertebra L1:20:48 \
  --output_path derivatives/fea --threads 4 --no-solve
```

Inspect `*_qc_3d.webp`: the body should be upright, processes posterior, and
both PMMA caps in contact with the body. Check `*_bc_audit.json` and
`*_qc_metrics.csv`. Registration/segmentation failures should be resolved
before interpreting mechanical results.

## Solve and report both regions

```bash
ogoFEA spine qct.nii.gz labels.nii.gz --vertebra L1:20:48 \
  --pistoia_mask_label 20 --critical_volume 12 --masked_critical_volume 35 \
  --critical_strain 0.007 --output_path derivatives/fea --threads 4
```

This reports reaction force at 0.68%, stiffness, full-bone Pistoia at 12%
critical volume, and body-only Pistoia at 35%, all in `*_results.csv`.
Critical volume is a percentage of evaluated bone; critical strain 0.007 means
0.7% effective strain. These criteria must be justified for the intended study.
The model-space ROI is saved as `*_pistoia_mask.nii.gz`.

For a separate binary body mask, replace `--pistoia_mask_label 20` with
`--pistoia_mask body_mask.nii.gz`. The full-bone result is still included.
For full-bone Pistoia alone, omit the ROI and add `--run_pistoia`.

## Review a cohort

```bash
ogoValidateFEA derivatives/fea --output derivatives/qc --gallery-zip --workers 4
```

Extract `gallery.zip`, open `gallery.html`, review and export the inclusion CSV.
Browser decisions do not overwrite project files. See the
[workflow guide](../../ogo/fea/README.md) for defaults, outputs and legacy models.
