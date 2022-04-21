# Usage

## Individual inclusion

Wrapper script [MCM_wrapper.sh](src/MCM_wrapper.sh):
* reads the [configuration file](config/config) 
* uses the files in [data/raw](data/raw) 
* includes the variants and individuals for this study 

The following scripts are invoked by the [wrapper](src/MCM_wrapper.sh):
* [prep_SNPs.sh](src/prep_SNPs.sh)
    Per cardiomyopathy (CM) this script takes the downloaded files of SNPs from ClinVar and VKGL and outputs a preliminary SNP list.
    Before executing this script, perform the following steps to obtain the input files:
    * In [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), search for all phenotype names
    * Filter for Pathogenic and Likely Pathogenic variants
    * Download the file with these features:
        * Format: Tabular (text)
        * Sort by: Location
    * Save the file (in this case: [data/raw/ACM_clinvar_result_LP.txt](data/raw/ACM_clinvar_result_LP.txt))
    * Download all [VKGL](https://vkgl.molgeniscloud.org/menu/main/dataexplorer?entity=vkgl_public_consensus&hideselect=true&mod=data#) results for the selected genes (in this case: [data/raw/ACM_VKGL.txt](data/raw/ACM_VKGL.txt))
* [extract_IID.sh](src/extract_IID.sh)
    * [Extract_IID_WES_UKB.R](src/Extract_IID_WES_UKB.R)
* [phenish.sh](src/phenish.sh)
* [Match_controls.R](src/Match_controls.R)
* [clean_pheno.sh](src/clean_pheno.sh)
    * [Combine_pheno.R](src/Combine_pheno.R)
    * [MCM_gene_summary.R](src/MCM_gene_summary.R)
    * [MCM_pheno_clean.R](src/MCM_pheno_clean.R)
The final dataset is then stored as data/processed/MCM_final_pheno.tsv and this is the input for all analyses and created tables and figures.
