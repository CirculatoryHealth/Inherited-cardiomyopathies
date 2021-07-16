# Inherited Cardiomyopathies
This project aimed to assess the prevalence of pathogenic mutations associated with ACM, DCM or HCM in the general UK population using UK Biobank data and assess the association of these mutations with cardiovascular related phenotypes.

## Curated genes
For this study, a curated list of genes was used to filter the (likely) pathogenic mutations.


## Pathogenic mutations of cardiomyopathies

### [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/)   
Here we extract [likely] pathogenic mutations for cardiomyopathies (ACM, DCM and HCM) in ClinVar.

- In ClinVar search for the cardiomyopathy of interest:  
  1. ACM: ""Arrhythmogenic right ventricular cardiomyopathy"[dis] OR "Arrhythmogenic right ventricular dysplasia"[dis] OR "Arrhythmogenic cardiomyopathy"[dis]"
  2. DCM: "Dilated cardiomyopathy"[dis]
  3. HCM: "Hypertrophic cardiomyopathy"[dis]

- Filter for pathogenic and likely pathogenic on the left side list.
- Download file with these features; Format:Tabular (text), Sort by:Location.

### [VKGL database](https://vkgl.molgeniscloud.org/menu/main/dataexplorer?entity=vkgl_public_consensus&hideselect=true&mod=data#)
The VKGL database was also used to extract (likely) pathogenic mutations for the cardiomyopathies. They were extracted using the following steps:
- In the data item filters wizard:
  - Select Classification: (Likely) Pathogenic
  - In the Gene box, insert the genes from the curated list per cardiomyopathy
  - Download file with these features: CSV


### Prepare the final list of (likely) pathogenic mutations
- Run the `prep_SNP.sh` script with the following arguments:
  - Argument 1: The phenotype for which the ClinVar-file was downloaded, so **ACM, DCM, or HCM**
  - Argument 2: The name of the ClinVar-file
  - Argument 3: The name of the VKGL-file
  - This script expects the ClinVar and VKGL in the `data/temp` folder
  - It also expects a file called `HCM/DCM/ACM_genes.txt` in this folder containing one gene per line of the curated gene list for this cardiomyopathy
```
prep_SNPs.sh [arg1] [arg2] [arg3]
```

- Include indel variants to extracted data by performing the following in 2 terminal screens:
  - In screen 1: find overlap of insertion and deletion:    
    1. `bin/overlap.pl ${DIS}_overlap_LP_WES_SNPs.txt 1 ${DIS}_overlap_LP_WES_SNPs.txt_position 1 -v | awk '$1 ~ /I/ {print $0}' | less` # get non-overlap Insertion    
    2. `bin/overlap.pl ${DIS}_overlap_LP_WES_SNPs.txt 1 ${DIS}_overlap_LP_WES_SNPs.txt_position 1 -v | awk '$1 ~ /D/ {print $0}' | less`  # get non-overlap Deletion   
  - In screen 2: search in `${DIS}_LP_positionID`:   
    1. Open the file and search each position from screen 1 "column2", then check if the indel is the same   
    2. Compare indel in overlap to the one in `${DIS}_LP_positionID` file
    3. Add matched indels to a file and save

:warning: FOR ACM, **Mimount** provides extra SNPs for both GRCh37 and GRCh38 in `WES_ACM_mB.csv`. We run the following command to obtain all SNPs in GRCh38 and add to the the list.
```
bash Mimount_ACM_SNPs_script.sh
```

## Extract UKB-participants carrying a mutation
After completion of the pathogenic mutation list, the individuals carrying a mutation are extracted using the following script and arguments"
- Argument 1: The phenotype, so **ACM, DCM, or HCM**
- Argument 2: The path to and name of the indel-file
```
extract_IID.sh [arg1] [arg2]
```


##
