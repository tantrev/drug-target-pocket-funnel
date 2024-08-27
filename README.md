# drug-target-pocket-funnel

This repository contains the scripts, data, and results from the paper "In Silico Target Identification Reveals IL12B as a High-Potential Candidate for Small Molecule Drug Development". It serves as a record of the computational methods used to identify IL12B as a potential drug target for small molecule development.

Contents:
- Custom scripts used for data processing and analysis
- Databases utilized in the filtering and analysis process
- Intermediate and final processed data
- Nextflow pipeline for specific parts of the analysis
- FTMap simulation results of the 1F42 PDB structure

## Technical Details

- The Nextflow script uses absolute paths and was run with Nextflow version NXF_VER=22.04.4
- Miniconda was used to install Python, dependent packages, and programs required for running the Nextflow pipeline

## Disclaimer

Complete reproduction of the paper's pipeline involves manual steps and requires the Schrödinger Small-Molecule Drug Discovery Suite (and an appropriate license). These are not included in this repository.

Please note that some scripts in this repository may not function properly without the appropriate "sitemap_results.csv", "sitemap_results_with_quality_metrics.csv", and "final_filtered_sites.csv" files. These files contain Schrödinger-specific outputs that have been excluded from the public repository to comply with the terms of Schrödinger's End User License Agreement. The affected scripts are provided for reference purposes only.

Accessing the complete data and code may require a Schrödinger software license. We have endeavored to provide all other non-Schrödinger files necessary to understand and reproduce our work.

Note: This repository is provided for transparency and reproducibility purposes related to the published paper. The scripts and data are specific to this project and its unique workflow.