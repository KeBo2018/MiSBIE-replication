# MiSBIE analysis data and code

Data and analysis scripts supporting *Mitochondrial Energy Transformation Capacity Influences Brain Activation During Sensory, Affective, and Cognitive Tasks*.

This package reproduces the principal statistical results from the non-image analysis data. It includes the final Figure 2 paired-response and Figure 3 scripts, portable Figure 1 and Figure 4 analyses, supplementary analyses, two-group power calculations, and historical imaging source code. Large participant imaging files are excluded.

## Run in MATLAB

1. Download the complete ZIP and use **Extract All**. Copying only the `.m` files does not include the inputs.
2. Inside the extracted `MiSBIE_data_code` folder, open **`START_HERE.m`** in MATLAB and click **Run**. It checks the inputs and then runs the analyses. To check the inputs alone, run `code/CHECK_MiSBIE_Setup.m`.
3. Keep `code/` and `data/` together. Each `RUN_MiSBIE_*.m` script can also be run separately. If a script has been moved, it opens a folder selector: choose the extracted `MiSBIE_data_code` folder (its `data` or `code` folder also works). This location is remembered for subsequent scripts.
4. Results are saved in the selected package's `outputs/` folder. No Excel selection or editing of local drive paths is required.

```text
MiSBIE_data_code/
    START_HERE.m
    code/
        RUN_ALL.m
        CHECK_MiSBIE_Setup.m
        RUN_MiSBIE_Figure1.m
        ...
    data/
        participants/analysis_data.csv
        phenotype/phenotype_scores.csv
        figure2/Paired_condition_scores.mat
        ...
```

**Previously copied scripts:** open the new files directly in MATLAB to avoid running an old copy on the MATLAB path. `which RUN_MiSBIE_Figure1 -all` lists competing copies. A path error mentioning `RevisionAfterPNAS/data` means the earlier script expected the released data there; it does not mean those files exist in the historical analysis folder. The path-fix release recognizes the complete extracted package or asks you to locate it. It cannot supply data from a ZIP that has not been extracted. See [docs/PATH_FIX.md](docs/PATH_FIX.md).

Use MATLAB R2020b or later with Statistics and Machine Learning Toolbox. Install CANlab Core and MediationToolbox to reproduce the original mediation inference. Without `mediation.m`, the Figure 3 script explicitly reports its percentile-bootstrap alternative; its significance annotations are not equivalent to CANlab's. CANlab is also needed for the optional Figure 1 permutation test. Dependencies and version limitations are listed in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md).

## Figure 1 and k-means repair (0.1.2)

The Figure 1 input CSVs contain the original 110 participants and 202 numeric variables, in the original order. All values and missingness match `Ke_DataforCorrelationMatrix_Cleaned_CK.xlsx`; only participant identifiers were replaced. The additional MRI metadata workbook is replaced by the bundled, ID-linked `cohort_groups.csv`. All figure inputs are already under `data/`.

All direct study k-means calls use the supplied **`kmeans_Matlab.m`**, included in `code/`, with the original options. Figure 1 itself uses Ward hierarchical clustering. The setup check exercises the uniquely named helper before running the figures. See [docs/FIGURE1_KMEANS_FIX.md](docs/FIGURE1_KMEANS_FIX.md).

## Contents

| Directory | Contents |
| --- | --- |
| `data/participants` | Release IDs, task-response scores, behavioral eligibility, age, NMDAS, GDF15, motion, demographic and derived structural measures |
| `data/phenotype` | 202 phenotype variables and 134 additional numeric clinical fields for the 110-person cohort |
| `data/figure2` | Paired condition scores, contrast scores across 100 saved CV repetitions, classification accuracy and effect sizes |
| `data/summary` | Tables S1–S5, full-precision Table S1 inputs, and archival Figure 1 effect-size selections |
| `data/figure_assets` | Existing Figure 2 map and schematic background; this artwork is not a new imaging analysis |
| `code` | Standalone MATLAB scripts and an independent Python numerical verification |
| `code/legacy_analysis` | Historical analysis scripts, with local home paths and literal study IDs redacted; reference only |
| `code/legacy_imaging` | Historical preprocessing, SVM fitting and resting-state scripts; original imaging inputs and path configuration required |
| `reference_results` | Independently calculated values and archived reference matrices for checking new runs |
| `docs` | Data dictionary, result-to-file map, provenance, validation and release notes |

## Independent numerical verification

From this directory, using Python 3.10 or later:

```sh
python -m pip install -r code/requirements.txt
python code/verify_results.py
```

The verification checks the task sample sizes, age-adjusted NMDAS correlations and Bayes factors, brain–behavior correlations, subgroup comparisons, OLS mediation coefficients, task effect sizes, GDF15 results, Table S1, the 304 correlations in each of Tables S2 and S3, and two-group power. It writes `outputs/verification/`.

MATLAB was not available during package preparation. The independent numerical verification passed; native MATLAB graphics and CANlab bootstrap inference have not been executed for this release. Run MATLAB before tagging a final public release.

## Reproduction boundaries

- The final Table S1 has 214 measures. Participant data validate 194; 20 retain original full-precision summary statistics because the corresponding raw variables were not found. This is explicit in the table's provenance column.
- All 304 coefficients, p values and sample sizes in each of Tables S2 and S3 reproduce to reporting precision. The historical `fdr()` output does not numerically equal standard Benjamini–Hochberg q values. Both are retained; the q < .05 classifications are unchanged. See [docs/VALIDATION.md](docs/VALIDATION.md).
- Several archived Figure 1A selection files exist. The exact final selection has not been uniquely established. The archived selections remain in data/summary. Figure 1 clustering now follows the user's supplied `Figure1_PartialCorr_SignAligned_k11_FINAL.m`: sign alignment, participant-row interpolation, partial Spearman correlation, Ward clustering at k=11, original labels/colors, and seeded t-SNE. Random t-SNE/permutation outputs can vary with software versions.
- Figure 2's small response plots recompute from saved held-out scores. The brain-map background is retained artwork. Re-estimating activation/SVM maps and resting-state connectivity requires excluded imaging inputs.
- NODDI sensitivity summaries were not located. Reported Tables S4–S5 and the historical source are included; the portable supplementary script computes the available covariate analyses and reports the missing NODDI input.

These boundaries are mapped by figure and table in [docs/RESULTS_MAP.md](docs/RESULTS_MAP.md). No direct between-task slope comparisons have been added.

## Participant data and public release

Original study identifiers were replaced with random release identifiers shared across the participant CSV files. Free-text clinical notes and event labels are excluded. Figure 2 `PairID` values are local to each task and must not be treated as links to participant `ReleaseID` values. The private ID key is not included in this repository or its ZIP.

The numeric data still describe a rare clinical cohort and include age, genotype categories and clinical measures. Removal of identifiers is not a determination that unrestricted public release is permitted. The study team must confirm the consent/institutional basis for public participant-level sharing before making the repository public. If public release is not permitted, retain these data in an approved controlled-access repository and adapt the data-availability statement.

Publication DOI and author-approved code/data licenses have not been assigned. See `LICENSE_NOTICE.txt` and [docs/UPLOAD_GUIDE.md](docs/UPLOAD_GUIDE.md). Participant identifiers and the private linkage key are not included in this repository.
