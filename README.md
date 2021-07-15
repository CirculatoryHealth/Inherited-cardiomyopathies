# Inherited-cardiomyopathies
This project aimed to to assess the prevalence of pathogenic mutations associated with ACM, DCM or HCM in the general UK population using UK Biobank data.


## Pathogenic mutations of cardiomyopathies in ClinVar  
Here we perform extraction of [Likely] pathogenic mutations for cardiomyopathies (ARVC, DCM and HCM) in ClinVar.
### [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/)   
- In ClinVar search for each cardiomyopathy:  
  1. ARVC: ""Arrhythmogenic right ventricular cardiomyopathy"[dis] OR "Arrhythmogenic right ventricular dysplasia"[dis] OR "Arrhythmogenic cardiomyopathy"[dis]"
  2. DCM: "Dilated cardiomyopathy"[dis]
  3. HCM: "Hypertrophic cardiomyopathy"[dis]

- Click on pathogenic and Likely pathogenic on the left side list.
- Download file with these features; Format:Tabular (text), Sort by:Location.
- Upload file to HPC and rename it with DCM/ARVC/HCM__P_and_LP_v3.txt.
- Run following command with [Argument1] which must be **ARVC, DCM, or HCM**:
```
prep_clinvar_SNPs_script1.sh [arg1]
```
- Include Indel variants to extracted data by performing the following in 2 terminal screens:
  - In screen 1: find overlap of insertion and deletion:    
    1. `/hpc/local/CentOS7/dhl_ec/software/overlap.pl ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt 1 ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /I/ {print $0}' | less` # get non-overlap Insertion    
    2. `/hpc/local/CentOS7/dhl_ec/software/overlap.pl ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt 1 ${TEMP}/${DIS}_overlap_path_SNPs_WES_SNPs.txt_position 1 -v | awk '$1 ~ /D/ {print $0}' | less`  # get non-overlap Deletion   
  - In screen 2: search in `${DIS}_P_and_LP_positionID`:   
    1. Use `less` to view the file then use "/" to paste each position from screen 1 "column2", then check if the indel is the same.   
    2. Compare indel in overlap to the one in `${DIS}_P_and_LP_positionID` file, then write them in excel sheet.
    3. Add matched indels to `overlap_path_SNPs_WES_SNPs.txt`.   
