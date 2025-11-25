# BinderBinderAnalyser
IgM binding array processing pipeline for the Binder lab from input raw csv files to quantification of feature importance

## Run on CeMM cluster
snakemake --profile ~/cemm.slurm.sm/

### Metadata description
1. Annotation table:
   List of samples (identical to sample names), comma-separated with columns providing information about the tissue and organism (not used in the pipeline, but relevant for potential downstream analysis)
2. Metadata table
Metadata table should contain fixed first six colummns with sample annotaiton introduction:
Sample_full_ID - sample id, same as in file names and anntation table, followed by columns providing information about the patient encrollment in the cohort: Record ID,Event Name,IRB Protocol,Date of Study Enrollment ,Age at Consent. This information is not used downstream, however those columns are skipped by the pipeline, so it's important to keep them (even as empty columns).
Further columns contain numeric and categorical metadata. If there is no information avaliable for a given sample, please keep the column empty! Should be comma-separated
3. Peptide annotation table
For each peptide coordinate prepare a tab-separated information file, contating the following information:
- row, column (see the input data description)
- Coordinate (row_column)
- sequence (unique aminoacid sequence of each pepetide)
- library (which peptide library is profiled - here we expect to have two peptide libraries)
- is_duplicate (yes/no column. It is recommended to include several peptides more than ones as control. If a peptide is added a second time, indicate it with a yes. However first instance of each peptide should be indicated as no)
- sequence_label (used for plotting)


### Input data description
CSV files named <sample>.csv are processed output of the ImageJ software.

Each CSV file contains two 

## Analysis steps
-Quality control and assesment of technical replicates
-Quantile normalisation
-Correlation of data with clinical parameters
-Random-forest-based evaluation of feature (peptide) importance