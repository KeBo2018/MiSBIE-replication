# Software dependencies

## Included-data MATLAB analyses

- MATLAB R2020b or later.
- Statistics and Machine Learning Toolbox for partial correlations, clustering, t-SNE, noncentral-t power, regression and density estimation.
- Figure 2 alone uses its own density helper and only base MATLAB graphics/numerical functions.
- [CANlab Core](https://github.com/canlab/CanlabCore) and [MediationToolbox](https://github.com/canlab/MediationToolbox) for the original `mediation` bootstrap inference and the optional `clusterdata_permtest` routine. Install these independently and follow their setup instructions. They are not bundled or relicensed here.

The user-supplied `code/kmeans_Matlab.m` is included unchanged (MathWorks copyright 1993–2020 retained). It still requires Statistics and Machine Learning Toolbox and its internal functions; the rename prevents accidental resolution to CANlab kmeans. It is third-party code and is not covered by any new license for the study code. Bioinformatics Toolbox is optional for the original Figure 1 extra `mafdr` diagnostic; when unavailable, its adjusted values are marked uncalculated and the main figures still run.

No MATLAB executable was available during preparation. The MATLAB scripts were adapted and inspected, but native execution is still required before a final version tag. The independently executed Python verification is not a replacement for that check.

The original CANlab/MATLAB commit or release versions were not recorded in the available source files. Use the authors' original environment if available and record `ver`, `which mediation -all`, `which clusterdata_permtest -all`, and toolbox Git commit IDs with the final rerun. Do not claim bitwise reproduction from an unpinned current toolbox release.

## Python verification

Python 3.10 or later with NumPy and SciPy; install from `code/requirements.txt`. No spreadsheet-reading package is needed because release inputs are CSV/MAT. The local verification environment is recorded in `reference_results/environment.json`.

## Historical imaging workflows

The archived scripts additionally use SPM, CANlab image objects, preprocessing outputs, masks/atlases, and in some cases signal-processing utilities. [SPM software](https://www.fil.ion.ucl.ac.uk/spm/software/) and the [CANlab dependency instructions](https://github.com/canlab/CanlabCore#dependencies-these-should-be-installed-to-use-this-toolbox) provide upstream installation information.

Subject imaging files, task events/logs, atlases and preprocessing derivatives are not included. The source paths indicate historical input filenames; they are not executable configurations for a new machine. Preserve the original preprocessing version and contrasts when rerunning them.
