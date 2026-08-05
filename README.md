# Single Nucleus RNA Sequencing (snRNAseq) Processing and Analysis

## Overview
This repository contains steps for processing and analalyzing single nucleus RNA sequencing (snRNAseg) data and can also be applied to single cell RNA sequencing (scRNAseq) data. Pre-processing steps are not included in this repo, but we recommend the use of CellRanger for preprocessing and initial QC, in adherence with best practices. The steps provided in this repo expect pre-processed sample Adata files (output from CellRanger) as input. Both source code (Python and R scripts) containing processing/analysis logic and SnakeMake files to run steps in workflow (individually using STANDALONE mode or as a full workflow) are provided.

## Workflow Description 
The workflow contains all steps for processing pre-processed snRNAseq data, generating analysis ready adata files, and running standard downstream analyses including expression analysis, differential expression analysis, and cell type enrichment analysis.

### Processing Steps
1. Merge sample adata files by region.
    - Inputs: Sample Adata Files (stored in specified directory)
    - Outputs: Region-stratified, merged Adata files containing all samples from specified region.

    1.1.  List sample adata files
        - Generates a table of sample IDs with associated Adata file path and sample region label

    1.2. - Merge sample adata files 
        - Genrates region-stratified, merged Adata files, using the sample Adata table generated in step 1.1.

2. QC filtering
    - Inputs: Region-stratified, merged Adata files
    - Outputs: Filtered region-stratified, merged Adata files
    - Default QC metrics:
        - N Genes by Count
        - Total Counts
        - % Counts Mitochondria
        - Doublet Score
    - Recommended QC thresholds:
        - N Genes by Count in range 300 - 10000
        - Total Counts in range 500 - 100000
        - % Counts Mitochondria less than 10
        - Doublet Score less than 0.2

    2.1. Plot QC metrics
        - Generates violin plots for QC metrics of interest that will be used for filtering
    
    2.2. Filter on QC Thresholds
        - Generates filtered region-stratified, merged Adata files based on specified QC metric thresholds

3. Metadata annotation
    - Inputs: Filtered region-stratified, merged Adata files
    - Outputs: Filtered region-stratified, merged Adata files with full sample metadata annotations

4. Cell type mapping w/ Allen's CellTypeMapper (Map My Cell backend)
    - Inputs:
        - Filtered region-stratified, merged Adata files with full sample metadata annotations
        - Reference taxonomy files: Precomputed stats and cell type marker genes
        - MMC gene mapper database (can be downloaded via ABCcache)
    - Outputs:
        - MMC extended results JSON per region
        - MMC results CSV per region
        - MMC log per region
        - Filtered region-stratified, merged Adata files with full sample metadata annotations and cell type mapping labels

    4.1 Run MMC
        - Generate MMC output files

    4.2 Add cell type labels
        - Generate Filtered region-stratified, merged Adata files with full sample metadata annotations and cell type mapping labels

5. Feature Selection
    - Inputs:
        - Filtered region-stratified, merged Adata files with full sample metadata annotations and cell type mapping labels
        - MMC marker genes file
        - [Optional] GWAS summary stats, specifying genes of interest
    - Outputs:
        - Reduced filtered region-stratified, merged Adata files with full sample metadata annotations and cell type mapping labels
        - All genes CSV
        - Reduced genes CSV
    - Notes: Retains highly variable genes (HVGs) + Cell Type Marker Genes + Genes of Interest (from GWAS summary statistics)

6. Integration 
    - Inputs:
        - Reduced filtered region-stratified, merged Adata files with full sample metadata annotations and cell type mapping labels
    - Outputs:
        - Reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, dimensionality reduction embeddings, and cell type prediction labels
    - Notes:
        - scVI trains a latent space model based on expression counts, effectively reducing dimensionality and corrected for batch effects
        - scANVI utilizes the scVI latent space model to train a cell type prediction model based on cell type mapping labels. Used to predict cell type labels for unresolved cells.

    5.1 scVI integration (dimensionality reduction)
        - Generates reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, and dimensionality reduction embeddings

    5.2 scANVI integration (cell type prediction)
        - Generates reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, dimensionality reduction embeddings, and cell type prediction labels

7. Clustering 
    - Inputs: Reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, dimensionality reduction embeddings, and cell type prediction labels
    - Outputs: Reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, dimensionality reduction embeddings, cell type prediction labels, and UMAP embeddings

8. Prepare final Adata files
    - Inputs: Reduced filtered region-stratified, merged Adata files with full sample metadata annotations, cell type mapping labels, dimensionality reduction embeddings, cell type prediction labels, and UMAP embeddings
    - Outputs: Final Adata files, comprising of previously specified attributes along with cell type assignment labels and log1p normalized expression counts (raw counts stored in "counts" layer)
    - Notes: Cell type assignment labels prioritize MMC cell type mapping labels for resolved cells and fall back on scANVI cell type prediction labels for unresolved cells.

9. Report generation
    - Inputs: Final Adata files, comprising of previously specified attributes along with cell type assignment labels and log1p normalized expression counts (raw counts stored in "counts" layer)
    - Outputs:
        - UMAP plots for specified groups (default groups: MMC cell type mapping labels, scANVI cell type prediction labels, cell type assignments, batch, leiden_res_0.05, leiden_res_0.10, leiden_res_0.20, leiden_res_0.40) and features (default features = default QC metrics)
        - Artifact metric plots and results CSV
    

### Downstream Analysis
1. Cell type expression analysis 
    - Inputs: Final adata object generated by processing workflow
    - Outputs: 
        - Region specific cell type x gene expression plots (dotplots)
        - Region specific cell type x gene "expressed" tables. Gene is expressed if...
            - Mean log1p norm expression for gene across all cells of a specific type is greater than 0.25
            - Fraction of cells of specific type showing non-zero expression for gene is greater than 0.1

2. Differential expression analysis (Case vs. Control)
    - Inputs: Final adata object generated by processing workflow
    - Outputs:
        - Region specific differential expression results
        - All significant differential expression results across regions

3. Cell type enrichment analysis
    - Inputs:
        - Final adata object generated by processing workflow
        - GWAS summary statistics OR precomputed MAGMA GWAS files
    - Outputs:
        - Region specific cell type enrichment results
        - Region specific cell type enrichment bar plots
    - Notes: GWAS results can provided as summary statistics, which will be transformed in MAGMA ready files (use workflow mode = "full"), or as precomputed MAGMA GWAS files (use workflow mode = "precomputed")

## Repository Structure
- base/
    - config.yml - used to define workflow parameters that are fed into scripts. Template contains default values, that can be updated by user.
    - Snakefile - SnakeFile used to run full workflow. Used STANDALONE=FALSE to run full workflow.
    - /env/ - environment files for Docker container. Can be used to build a conda env.
    - /src/ - source code scripts.
    - /modules/ - SnakeFiles for each workflow step. Use STANDALONE=TRUE to run individual step modules.

## Environment
Provided in this repo are scripts to build a reproducible environment containing all necessary packages and libraries to run each step in the processing/analysis workflow along with a Dockerfile to build an image for this environment and run the workflow within a container. To build a conda environment containing all necessary packages and libraries, follow these steps (ensure conda is installed):

1. conda env create -f env/environment.yml
    - Builds the conda environment and installs necessary packages (using conda and pip)
2. conda activate snRNAseq
    - Activates the environment
3. Rscript install_libraries.r
    - Installs necessary R libraries within the conda environment
4. git clone https://github.com/AllenInstitute/cell_type_mapper.git \
    && cd cell_type_mapper \
    && pip install .
    - Installs cell_type_mapper from source
5. snakemake -s Snakefile -j1 -p
    - Runs full processing and analysis workflow.
6. snakemake -s modules/<module.smk> -j1 -p
    - Runs individual specified step module.

To build the Docker image and run this workflow (or individual steps) within a container, follow these steps (ensure Docker is installed):

1. docker build Dockerfile .
    - Builds the Docker image containing all necessary libraries and packages.
2. docker run --rm \
  -v </path/to/data>:/snRNAseq/data \
  -v </path/to/results>:/snRNAseq/results \
  snrnaseq-base --cores <all>
    - Runs full processing and analysis workflow.
    - Make sure to specifify the paths to your data and results directories, along with the number of cores.
3. docker run --rm \
  -v </path/to/data>:/snRNAseq/data \
  -v </path/to/results>:/snRNAseq/results \
  snrnaseq-base \
  --snakefile modules/<module.smk> --cores <all>
  - Runs individual step module.
  - Make sure to specifify the paths to your data and results directories, along with the number of cores and the filename for the module you wish to run.