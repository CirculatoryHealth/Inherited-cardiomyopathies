# Usage

## Individual inclusion

Wrapper script [MCM_wrapper.sh](MCM_wrapper.sh):
* reads the [configuration file](../config/config) 
* uses the files in [data/raw](../data/raw) 
* includes the variants and individuals for this study 

The following scripts are invoked by the [wrapper](MCM_wrapper.sh):
* [prep_SNPs.sh](prep_SNPs.sh)
  Per cardiomyopathy (CM) this script takes the downloaded files of SNPs from ClinVar and VKGL and outputs a preliminary SNP list.
    Before executing this script, perform the following steps to obtain the input files:
    * In [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), search for all phenotype names
    * Filter for Pathogenic and Likely Pathogenic variants
    * Download the file with these features:
        * Format: Tabular (text)
        * Sort by: Location
    * Save the file (in this case: [ACM_clinvar_result_LP.txt](../data/raw/ACM_clinvar_result_LP.txt))
    * Download all [VKGL](https://vkgl.molgeniscloud.org/menu/main/dataexplorer?entity=vkgl_public_consensus&hideselect=true&mod=data#) results for the selected genes (in this case: [ACM_VKGL.txt](../data/raw/ACM_VKGL.txt))
* [extract_IID.sh](extract_IID.sh)
  Before executing this script, indels should be manually searched and included. In this case, this file is saved as [ACM_indels.txt](../data/raw/ACM_indels.txt). Also, SNPs from the VKGL database were checked in ClinVar to see if there was evidence for this SNP being associated with (one of the) cardiomyopathies. If there was no evidence, the SNP was added to [ACM_remove_overlap.txt](../data/raw/ACM_remove_overlap.txt). This script finalizes SNP inclusion and extracts the UKB participants that carry these variants. It invokes the following R-script to finalize the individual list per CM. 
    * [Extract_IID_WES_UKB.R](Extract_IID_WES_UKB.R)
      This script makes a column per gene showing the carrier status for each individual.
* [phenish.sh](phenish.sh)
  This script takes, per CM, the output of the previous script, a fam-file of the UKB WES data and the complete CMR dataset and extracts the required phenotypes from the UKB, merging it with the CMR data.
* [Match_controls.R](Match_controls.R)
  Taking the output files of all CMs of the previous script, it matches controls to the selected cases. The potential control population is the full WES population, excluding the cases, and matching is performed based on age, sex, ethnicity and CMR availability. 
* [clean_pheno.sh](clean_pheno.sh)
  This script finalizes the phenotype file, combining all individuals, genotype and phenotype information. It does so, invoking the following scripts:
    * [Combine_pheno.R](Combine_pheno.R)
      This script takes the phenotype files of all groups, harmonizes the headers and combines them.
    * [MCM_gene_summary.R](MCM_gene_summary.R)
      This script summarizes the genotype information in few columns.
    * [MCM_pheno_clean.R](MCM_pheno_clean.R)
      Using all output of previous scripts, the outcomes are defined, dates of diagnoses are included, columns are renamed and some variables are organized.

  
The final dataset is then stored as MCM_final_pheno.tsv and this is the input for all analyses and created tables and figures.

# Analyses

The prevalence and its confidence interval were (manually) calculated by [CI_prevalence.py](CI_prevalence.py).
Difference testing is performed by the scripts [Create_TableS678.R](Create_TableS678.R) and [Create_TableS9.R](Create_TableS9.R). These files take the final dataset as first argument and the output file names as other arguments. These results are visualized by [Create_Figure4_S2_S3.R](Create_Figure4_S2_S3.R).
