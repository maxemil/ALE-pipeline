#! /bin/bash

display_usage() {
  echo "wrong number of arguments, usage:"
  echo "bash aggregat_node_tables.sh <file-nodelist> <file-events> <path-to-faa> <glob-to-eggnog-annotations> <prefix-species-tree> <path-output-tables> <file-tree-certainty>"
}

if [  $# -ne 7 ]
then
		display_usage
		exit 1
fi


STR=""
while IFS=$',' read -r -a myArray
do
        STR=$STR" ${myArray[1]} ${myArray[2]}"
done < $1

echo "anapy3 /local/two/Software/ALE-pipeline/bin/match_cluster_eggnog.py -f $3 -a $4 -s $5 -c cognames2003-2014.tab -e $2 -np $STR -t $7"
#
# mkdir -p $6/tab_files
# mkdir -p $6/xls_files
#
# while IFS=$',' read -r -a myArray
# do
#     ssconvert -M "${myArray[0]}.xls" "${myArray[2]}_cluster_OG.tab" "gains_${myArray[1]}to${myArray[2]}.tab" "losses_${myArray[1]}to${myArray[2]}.tab"
# done < $1
#
# mv *cluster_OG.tab losses*.tab gains*.tab $6/tab_files
# mv *.xls $6/xls_files
