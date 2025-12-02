# IgM assay analyser
IgM binding array processing pipeline for the Binder lab from input raw csv files to quantification of feature importance and code for the analysis, presented in the paper Deroissart et al., XXXX

## IgM processing pipeline

### Running the pipeline
To perform a test run from an environment with installed snakemake (tested on snakemake8).
```
snakemake --profile <your_HPC_profile.sm>
```
### Setup
To prepare the run you need to:  
- specify the project config path adjust the cluster settings to your HPC in the ``config/pipeline_config.yml``
- specify the paths and the run-specific parameters in the project config file. The example config file with explanations is provided in test.
The pipeline can run in two modes: either with or without preprocessing (from individual samples with technical replicas to QN files), regulated by flag skip_preprocessing. Flag should be set to True if only goal is to repeat downstream analysis with different parameters or reproduce downstream analysis on two cohorts merged together.
To provide a custom (i.e. batch-corrected from two runs) QN files use a flag use_custom_QN and in this case provide path to the QN file to use in custom_QN_file file.
Default pipeline separates samples by sex for downstream analysis. Make sure that the metadata file contains the *Sex* column, that specifies Male and Female for patient annotation. If other value groups are used, please modify it at *line 117 of the snakemake file*.

### Metadata description
1. Annotation table:
   List of samples (identical to sample names), comma-separated with columns providing information about the tissue and organism (not used in the pipeline, but relevant for potential downstream analysis). Example annotation table is provided in ``test/config``
2. Metadata table
Metadata table should contain fixed first six colummns with sample annotaiton introduction:
Sample_full_ID - sample id, same as in file names and anntation table, followed by columns providing information about the patient encrollment in the cohort: Record ID,Event Name,IRB Protocol,Date of Study Enrollment ,Age at Consent. This information is not used downstream, however those columns are skipped by the pipeline, so it's important to keep them (even as empty columns).
Further columns contain numeric and categorical metadata. If there is no information avaliable for a given sample, please keep the column empty. Make sure that column names are R-friendly and don't contain any special characters. 
File should be comma-separated. Example is provided in ``test/meta``
4. Peptide annotation table
For each peptide coordinate prepare a tab-separated information file, contating the following information:
- row, column (see the input data description)
- Coordinate (row_column)
- sequence (unique aminoacid sequence of each pepetide)
- library (which peptide library is profiled - here we expect to have two peptide libraries)
- is_duplicate (yes/no column. It is recommended to include several peptides more than ones as control. If a peptide is added a second time, indicate it with a yes. However first instance of each peptide should be indicated as no)
- sequence_label (used for plotting)
Example is provided in ``test/meta``.

### Input data description
CSV files named <sample>.csv are processed output from the ImageJ software.
Each CSV file contains two technical replicas with rows and columns corresponding to the coordinates of the plate in the assay. Example data is provided in ``test/data``

### Analysis steps
- Quality control and assesment of technical replicates
- Quantile normalisation
- Quantification of relationships between the data with clinical parameters (correlations, confounding analysis)
- Differential analysis based on the metadata column

## Downstream data analysis
In the paper 3 separate cohorts where used: 
CAVA cohort, profiled in 2 batches and BIKE cohorts for plasma and aorta plaque samples (matched). The two batches of CAVA cohort where processed independantly, then quantile normalised files with corrected for batch effect as shown in the ``downstream_analysis/XX.ipynb``. Created file was further re-analsyed within the pipeline, using the *custom_QN_file* option.
