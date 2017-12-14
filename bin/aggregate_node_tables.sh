#! /usr/bin/bash

STR=""
while IFS=$',' read -r -a myArray
do
        STR=$STR" ${myArray[1]} ${myArray[2]}"
done < $1

echo "bin/match_cluster_eggnog.py -f clusters/*.faa -a Emapper_annotations.tsv --species_tree species_tree_name -e events.txt -np $STR"

mkdir tab_files
mkdir xls_files

while IFS=$',' read -r -a myArray
do
    ssconvert -M "${myArray[0]}.xls" "${myArray[2]}_cluster_OG.tab" "gains_${myArray[1]}to${myArray[2]}.tab" "losses_${myArray[1]}to${myArray[2]}.tab"
done < $1

mv *cluster_OG.tab losses*.tab gains*.tab tab_files
mv *.xls xls_files
