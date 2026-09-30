# Upload and release

Upload only the contents of **MiSBIE_data_code**. Its sibling `private_preparation_do_not_upload` contains the confidential identifier lookup and local preparation audit and must stay outside GitHub. The supplied ZIP contains only MiSBIE_data_code.

1. Review public-sharing permission for the participant-level clinical data. The files are pseudonymized, not certified anonymous. If a controlled-access route is required, change the release contents and data-availability statement before publication.
2. Run `code/RUN_ALL.m` in MATLAB with the original toolboxes if available. Compare outputs against `reference_results` and review `docs/VALIDATION.md`.
3. Choose author/institution-approved licenses for the study code and data. No license has been assigned on the authors' behalf. Preserve third-party notices in historical files.
4. Create the GitHub repository and upload these contents, preserving `data/` and `code/`. Do not upload source clinical Excel workbooks, the parent sharing directory, or participant identifiers.
5. Add the repository URL, full author citation, publication DOI when available, and toolbox versions. Make a tagged release corresponding to the submitted/published analysis. A DOI can be added if the authors archive that release in a suitable repository.
6. Update the manuscript's data/code availability statement to describe the actual public and controlled-access contents. A generic GitHub profile is not the analysis repository URL.

`MANIFEST.csv` records file sizes and SHA256 hashes for the prepared package. `outputs/` is ignored by Git because it contains regenerable results. All source input values and important reference outputs are elsewhere in the repository.

This repository contains only the prepared release files. Original study files, the private identifier lookup, and local preparation audits are kept separately.
