#!/bin/bash

echo "Creating two similar files"
sed -i "" 's/ /_/g' *read.txt
head -1 CTV3.read.txt > RCT3.read.txt
awk 'NR > 1 {print $1, "NA", $2, "NA", $3, $4, $5, $6, $7}' RCT.read.txt | sed 's/ /\t/g' >> RCT3.read.txt

echo "Searching terms now..."
echo -e "$(head -1 CTV3.read.txt)\tFile" > GP_keyword_search.tsv

while IFS= read -r line; do

  echo "Searching for the following term now: ${line}"
  echo -e "$(grep -i ${line} RCT3.read.txt)\tRCT" >> GP_keyword_search.tsv
  echo -e "$(grep -i ${line} CTV3.read.txt)\tCTV" >> GP_keyword_search.tsv
  echo ""

done < "GP_verlanglijstje.txt"
