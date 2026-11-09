# X-Atlas/Orion derived panels and models

License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0),
https://creativecommons.org/licenses/by-nc-sa/4.0/

The files in this folder are derived from the X-Atlas/Orion genome-wide Perturb-seq atlas
(HCT116 and HEK293T cells; Huang et al., 2025), which is distributed under CC BY-NC-SA 4.0.
They are therefore released under the same license and kept apart from the K562-derived
deliverables in `deliverables/`, which carry the license of the main repository.
Non-commercial use only. Derivatives must be shared under the same license and credit both
X-Atlas/Orion and this work.

## Contents

- `panels/hct116_500_standard.tsv`, `panels/hek293t_500_standard.tsv`: greedy minimum mean squared
  error panels of 500 genes selected on the HCT116 (17,762 perturbations x 11,948 genes) and HEK293T
  (17,999 perturbations x 12,288 genes) atlases.
- `panels/pan_atlas_500_standard.tsv`: panel selected on the pooled reference below.
- `models/hct116_500_standard.bin`, `models/hek293t_500_standard.bin`: `perturb-impute` models
  fitted on each atlas (INFO and calibration from the grouped 10% inner holdout).
- `models/pan_atlas_500_standard.bin`: model fitted on K562 + HCT116 + HEK293T perturbations on
  the 6,414 genes measured in all atlases and in HepG2 cells. Rows of each atlas are scaled so
  that each atlas contributes equal weight to the reference covariance. This pooled reference
  recovered more of the specific response than any single atlas in RPE1, HepG2 and Jurkat cells.

Use with `perturb-impute impute --model <file> --input <panel z-scores> --out-dir <dir>`
(see the main README).
