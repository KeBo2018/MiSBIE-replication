# Folder lookup repair — 0.1.1-pathfix

The initial scripts assumed that the directory above the script contained `data/`. The reported errors searched `E:/Mito_DICOM/RevisionAfterPNAS/data`, where the packaged CSV/MAT files were absent. Those errors occur at input loading, before a statistical analysis is performed.

## Use the complete corrected package

1. Extract **all** files from `MiSBIE_data_code_GitHub_package_PathFix.zip` to any location, including the E: drive if desired.
2. Open `MiSBIE_data_code/START_HERE.m` from the extracted files and click Run. This opens the new entry point directly, avoiding a same-named old script on the MATLAB path.
3. The run starts with a setup check that finds the 12 required input files and loads the participant tables and Figure 2 MAT file. For a setup check alone, open `code/CHECK_MiSBIE_Setup.m`.

Each analysis also contains its own package locator. It checks the script location, nearby package folders, the current working directory, and a previously selected package location. If necessary, a folder selector asks for the extracted package. The root, `code`, or `data` folder can be selected. The chosen location is stored in the MATLAB preference `MiSBIE_Reproducibility/PackageRoot`. A missing or moved saved folder is rechecked rather than blindly reused.

Canceling selection stops with an explanation. An incomplete package stops with a list of missing files. Outputs are created under the selected package only after its required inputs have been located. No historical Excel workbook is substituted, and statistical routines, participant flags, source data and effect estimates are unchanged by this repair.

## Diagnose an old copy

```matlab
which RUN_MiSBIE_Figure1 -all
which RUN_MiSBIE_Figure2_Paired -all
which RUN_MiSBIE_Figure3 -all
which RUN_MiSBIE_Figure4 -all
```

If these show old copies, open the new script by its full path rather than invoking a same-named old copy at the Command Window. To forget a previously selected package:

```matlab
if ispref('MiSBIE_Reproducibility','PackageRoot')
    rmpref('MiSBIE_Reproducibility','PackageRoot');
end
```

The folder repair has been inspected and the numerical verification is rerun against a relocated extracted package. Native MATLAB execution remains unverified because MATLAB is unavailable on the preparation computer. This update addresses the reported missing-file failures; it does not claim that MATLAB graphics or toolbox-dependent analyses have been executed here.
