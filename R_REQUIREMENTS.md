# R and software requirements

The manuscript-generating session logs record the following environment.

## Tested environment

- R version: 4.5.0 (2025-04-11 ucrt)
- Platform: x86_64-w64-mingw32/x64
- Operating system: Windows 10 x64
- Locale: Chinese (Simplified) China UTF-8
- Time zone: Europe/Dublin

## Required R packages

Core analysis scripts require:

- `data.table`
- `ggplot2`

Recorded package versions in Step 5/6 session logs:

- `data.table` 1.17.0
- `ggplot2` 3.5.2

## Optional R packages

The manuscript table export script can use one of:

- `openxlsx`
- `writexl`

If neither package is installed, the final table script writes CSV outputs only and logs that XLSX export was skipped.

## PowerShell

Step 1 uses a Windows PowerShell script:

- `code/scripts/build_pcr_manifest.ps1`

The R scripts are written with project-root arguments and should be portable across operating systems after Step 1 outputs are available.

