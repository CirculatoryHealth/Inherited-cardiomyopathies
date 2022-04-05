#!/bin/bash
#
#SBATCH --job-name ukb_pheno                                     # the name of this script
#SBATCH --time=02:30:00                                          # --time=[max time, e.g. 02:02:01] - this is the time you think the script will take
#SBATCH --mem=1G                                                 # --mem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use

# Step 1: Download the DataDictionary from the UKBiobank site. (http://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=AccessingData)
# Step 2: Select the FieldID and Fiels you need for your research and copy them to your input file. Please take in mind that the header should be like the header in the exapmle below.
# Step 3: Make sure the space between the two columns is tab space as shown in the example below.
# Step 4: When you run the script, type the name of your inputfile after the command you use to run the script. For example: $ csh ukb_pheno.csh tab_file fieldID_file output_directory <project_name>

#							    FieldID Field
#							    3       Verbal interview duration
#							    4       Biometrics duration
#							    5       Sample collection duration
#							    6       Conclusion duration
#							    19      Heel ultrasound method

INPUT_TAB_FILE=$1
INPUT_FIELD_ID=$2
OUTPUT=$3
OUTPUT_FILE_NAME=$4

### Eternal Tools
SCRIPTDIR=/hpc/dhl_ec/mvanvugt/scripts
MERGE=/hpc/dhl_ec/aalasiri/ukb_pheno/merge_tables.pl
TRANSPOSE=${SCRIPTDIR}/transpose1.pl

# transpose header of the origin tab file
head -1 ${INPUT_TAB_FILE} > ${OUTPUT}/header_1_original.tab
${TRANSPOSE} ${OUTPUT}/header_1_original.tab > ${OUTPUT}/header_1_original.tab_trans

# column location for phenotypes
echo "Column.No FieldID" > ${OUTPUT}/pheno_location_v1
awk '{print $1,$1}' ${OUTPUT}/header_1_original.tab_trans | sed 's/\./ /' | sed 's/\./ /' | cat -n| awk '{print $1,$3}' | sort -k2 -u | sed '$d' >> ${OUTPUT}/pheno_location_v1
# number of measures
echo "NO.measure FieldID" > ${OUTPUT}/pheno_NO.measure_v1
awk '{print $1,$1}' ${OUTPUT}/header_1_original.tab_trans | sed 's/\./ /' | sed 's/\./ /' | cut -d' ' -f2 | uniq -c | awk '{print $1, $2}' >> ${OUTPUT}/pheno_NO.measure_v1

# phenotypes information file
echo "Extracting required phenotypes from ${INPUT_TAB_FILE}"
sed 's/ /\_/g' ${INPUT_FIELD_ID} > ${OUTPUT}/f.eid_modif_v1 #replace spaces between wards
# merge number of measures and input file
${MERGE} --file1 ${OUTPUT}/pheno_NO.measure_v1 --file2 ${OUTPUT}/f.eid_modif_v1 --index FieldID > ${OUTPUT}/output2.1_v1
# merge column location with previous
${MERGE} --file1 ${OUTPUT}/pheno_location_v1 --file2 ${OUTPUT}/output2.1_v1 --index FieldID | sed  's/NA/0/g' | awk '$3 != 0 && $4 != 0 {print $0}' | sort -k 1n > ${OUTPUT}/output2.2_v1

## Extract interest phenotypes from origin tab file
VAR=$(tail -n +2 ${OUTPUT}/output2.2_v1 | awk '{print $4,$4+$3-1}' | sed 's/ /-/g' | ${TRANSPOSE} | sed 's/\t/,/g')
cut -f1,$VAR ${INPUT_TAB_FILE} > ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.tab
echo "Just created a tab file containing your phenotypes in UKB: ${OUTPUT_FILE_NAME}_ukb_phenotypes.tab"

## Print output information file
echo "Processing the phenotypes information file"
echo "## Available phenotypes on ${INPUT_TAB_FILE}" > ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo  >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
cat ${OUTPUT}/output2.2_v1 >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info

echo  >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo  >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo  >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo "## Unavailable phenotypes  on ${INPUT_TAB_FILE}" >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo  >> ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info

awk '$3 == "NA" {print $1,$2}' ${OUTPUT}/output2.1_v1 >>  ${OUTPUT}/${OUTPUT_FILE_NAME}_ukb_phenotypes.info
echo "Just created an info file containing more details: ${OUTPUT_FILE_NAME}_ukb_phenotypes.info"

#####

rm ${OUTPUT}/header_1_original.tab
rm ${OUTPUT}/header_1_original.tab_trans
rm ${OUTPUT}/pheno_location_v1
rm ${OUTPUT}/pheno_NO.measure_v1
rm ${OUTPUT}/f.eid_modif_v1
rm ${OUTPUT}/output2.1_v1
rm ${OUTPUT}/output2.2_v1
echo "END"
