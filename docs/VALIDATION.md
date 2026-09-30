# Numerical validation and remaining boundaries

Prepared 30 September 2026 from the final manuscript analysis inputs. The original workbooks, code and manuscripts were preserved. Validation uses the exported release CSVs, not a hidden local workbook.

## Independently reproduced

| Check | Result |
| --- | --- |
| Phenotyping cohort | 110 unique release IDs; 40 patients and 70 controls |
| Primary N-back sample | 74 participants; 25 patients, 49 controls |
| Multisensory and cold-pain samples | 90 and 91 participants; 28 and 29 patients |
| Age-adjusted NMDAS/brain rho | N-back −0.6703603; multisensory −0.3705955; cold pain −0.1915210 |
| Within-task BF10, prior scale r=1 | 89.541921; 0.883389; 0.230625 |
| Mediation complete cases | 72; OLS path coefficients reproduced |
| Final-split task dz | 1.6787253; 2.1835211; 0.9648906 |
| Cold-pain display filter | Three pairs flagged; all 91 pairs retained in effect size and data |
| Figure 4 combined brain correlation | Group-adjusted Pearson r=−0.4847133, n=67 |
| Figure 4 combined behavior correlation | Group-adjusted Pearson r=−0.3281309, n=67 |
| GDF15 brain regression | df=63; GDF15 t=−2.81923, age t=−2.75149, group t=−2.95286 |
| Table S1 | 194 raw-variable matches, 20 inherited summaries; all 214 reconciled rows reproduced, 211-test BH family |
| Tables S2 and S3 | All 304 correlations, p values and sample sizes in each table match the reported precision |
| Two-group power | Recomputed with actual unequal group sizes and balanced required N per group |

`reference_results/` contains the numerical CSV outputs and machine-readable checks. `code/verify_results.py` regenerates them under `outputs/verification/`. Raw Spearman p values in the Python check use the asymptotic t formula; native MATLAB small-sample handling may differ. Partial correlations and deterministic coefficients are checked independently.

## Historical FDR values in Tables S2 and S3

The historical analysis calls an external `fdr()` function. Its saved `p_FDR` values differ from conventional monotone Benjamini–Hochberg q values for the same 304 p values, including entries saved as 1. The exact historical `fdr()` implementation/version was not found in the supplied code. The reported values are preserved in `data/summary`; the portable code explicitly exports both `ReportedFDRValue` and `QBH304Tests`.

This difference does not change any q < .05 classifications: Table S2 has 18 significant associations under either set of values; Table S3 has none. It does prevent a claim that the historical adjusted numerical values have been exactly regenerated. Resolve the function/version or approve a documented table update before claiming exact reproduction of those columns. No manuscript or table values were changed during packaging.

## What was not executed or established

- Native MATLAB execution, MATLAB plot layout and CANlab bootstrap inference. Inference from the original CANlab routine must not be conflated with the script's labeled fallback.
- Imaging preprocessing, voxelwise activation maps, SVM refitting and resting-state analyses. Large source imaging inputs are excluded at the authors' request; historical code and saved task-score outputs are included.
- Exact final Figure 1A variable selection. Several source versions remain; the archival wedge inputs are retained separately from the restored clustering script.
- Original Figure 1 random t-SNE/permutation realizations and toolbox versions. The supplied source uses participant-row interpolation; row order is preserved for that reason.
- Participant-level revalidation of the 20 missing Table S1 variables, or the unavailable NODDI sensitivity summaries.
- Ethical/institutional authorization for public release or an author-approved license. The local package has not been published.

## Small reporting distinction

The original Figure 4 plotting code specifies Spearman for the GDF15/NMDAS panel, unlike its Pearson brain/behavior panels. The available pre-exclusion N-back data yield rho=.75816 (22 patients), while the manuscript reports .75. The portable script exports Spearman and Pearson separately and does not change the manuscript.

## Figure 1 source confirmation and k-means repair

The user supplied `Figure1_PartialCorr_SignAligned_k11_FINAL.m` and `kmeans_Matlab.m`. The previous runnable Figure 1 adaptation omitted sign alignment; this was a substantive mismatch and is corrected in version 0.1.2. The original sign-aligned pipeline, k=11, palette and cluster labels are restored.

`figure1_source_check.json` verifies all 110 rows and 202 variables against the original workbook, including missingness and participant order. An independent numerical implementation identifies 107 flipped variables and cluster sizes 6, 8, 11, 14, 14, 15, 16, 18, 18, 28, 54, matching the multiset recorded in the supplied source. This does not assert identical MATLAB cluster numbering or t-SNE coordinates.

The k-means helper matches the supplied file byte-for-byte. All direct calls in the released study code use `kmeans_Matlab`; historical scripts already using that name retain it. The setup includes a fixed-start clustering check for execution in MATLAB. Native MATLAB/helper execution is still unverified here.
