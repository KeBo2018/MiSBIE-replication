# Figure 1 restoration and k-means name conflict

## Run

Extract the **complete** `MiSBIE_data_code_Figure1_KmeansFix.zip`. Open `code/RUN_MiSBIE_Figure1.m` from this new folder and click Run. Its inputs are already included:

- `data/phenotype/phenotype_scores.csv`: all 110 original participant rows and 202 numeric variables from `Ke_DataforCorrelationMatrix_Cleaned_CK.xlsx`.
- `data/participants/cohort_groups.csv`: ID-linked group labels for 40 patients and 70 controls. This replaces the original script's external metadata workbook.
- `code/kmeans_Matlab.m`: the exact supplied implementation, needed by the patient subgroup analyses, not by Figure 1's Ward clustering.

`code/CHECK_MiSBIE_Setup.m` checks the inputs and tests the supplied k-means helper. `START_HERE.m` runs the setup and all analyses. Keep `data/` and `code/` together. The folder selector/remembered package location remains available if a script has been moved.

## Restored Figure 1 workflow

The supplied `Figure1_PartialCorr_SignAligned_k11_FINAL.m` is the source for the runnable Figure 1 script. It retains numeric-variable selection, median/IQR standardization, sign alignment by the patient-minus-control mean, participant-row interpolation, group-adjusted partial Spearman correlations, Ward linkage and k=11. It also retains the original cluster names, palette, t-SNE options and seed, centroid labels, correlation heatmap and diagnostic figures.

The previous adaptation omitted sign alignment. With distance `sqrt(2*(1-r))`, changing correlation signs can change Ward clustering. The restored step is necessary for this replication. The source comment claiming that sign alignment cannot change clustering has been corrected; no additional clinical meaning is assigned to its data-derived sign convention.

All numeric values, missingness, variable order and participant row order match the source workbook. Only the participant identifier strings differ. Independent calculation gives 107 sign inversions and the same cluster-size multiset recorded in the source: 6, 8, 11, 14, 14, 15, 16, 18, 18, 28 and 54. Native MATLAB cluster numbering and t-SNE coordinates remain to be checked.

Outputs are saved under `outputs/Figure1/`, including the original panel names, PNG/PDF/FIG files, `partial_correlation_matrix.csv`, `phenotype_clusters.csv`, `sign_alignment_audit.csv`, participant order and the intermediate MAT file. The attached source covers the clustering panels; it does not establish the final Figure 1A effect-size selection. Archived Figure 1A inputs remain in `data/summary/`.

## Runtime corrections

- All input and output paths are relative to the selected package.
- The supplied `kmeans_Matlab.m` is copied unchanged. Direct study calls to `kmeans`/`k_means` use this unique name and the original options; scripts add their package code folder to the MATLAB path. Historical calls already using that name remain unchanged. CANlab itself is not edited.
- The original four-component text background color is replaced by a supported RGB white background. Cluster-border loops explicitly iterate over individual boundaries. The per-variable sign diagnostic prints each variable's own difference.
- The distance radicand is clamped at zero for roundoff at correlation 1.
- The optional CANlab permutation diagnostic retains k=2:40 and 100 permutations. It is seeded for future reproducibility and explicitly skipped if the dependency is absent. This does not fabricate an original unseeded permutation result.
- The extra NMDAS diagnostic retains the original `mafdr` call when available. Its default is Storey's procedure, although the source comment called it BH; the annotation is corrected. If the dependency is absent, raw correlations/p values are exported and the adjusted column is marked uncalculated. The manuscript's Tables S2/S3 are not changed by this diagnostic.

Relevant MATLAB behavior is documented in [fillmissing](https://www.mathworks.com/help/matlab/ref/fillmissing.html) and [mafdr](https://www.mathworks.com/help/bioinfo/ref/mafdr.html).

The supplied k-means implementation retains its MathWorks copyright notice and still requires Statistics and Machine Learning Toolbox. No new license for that third-party file is granted by this package. Native MATLAB execution was unavailable on the preparation computer; a successful data audit and Python calculation are not a claim that MATLAB graphics or toolbox execution have been tested.
