#!/bin/bash

YEAST_CYCLES=(cycle_alpha7 cycle_alpha42 cycle_alpha49)

for YEAST_CYCLE in ${YEAST_CYCLES[@]}; do
	IN="yeast_out/strand1/${YEAST_CYCLE}"
	OUT="yeast_abc4pwm_out/strand1/${YEAST_CYCLE}"

	for MLP in `ls ${IN}`; do
		echo " "
		echo " "
		echo " ########running abc4pwm on = $MLP"

		abc4pwm searching --pwm ${IN}/${MLP} --db_path SGD_yeast_cellCycle/ --output_directory ${OUT}/${MLP} --db_type folder --input_count 1 --db_count 0 --db_file_type '.mlp' --input_file_type '.mlp' --input_prob 0 --db_prob 1
		echo " ########DONE for $MLP"
	done

echo " "
echo "Combining scores for predicted mols"
FNAME_FILES=()
	for MY_FILE in `ls ${OUT} |grep mlp$`; do
		echo " ########running abc4pwm on = $MY_FILE"
		REQ_FILE="${OUT}/${MY_FILE}/search_results.txt"
		OUT_FILE="${REQ_FILE}.fname"
		sed  "s/$/\t$MY_FILE/" $REQ_FILE >$OUT_FILE
		echo " ########Written file: ${OUT_FILE}"
		FNAME_FILES+="$OUT_FILE "
	done
cat ${FNAME_FILES[@]} > "${OUT}/${YEAST_CYCLE}_matchted_mlps.txt"
sort -h -r -k 2 "${OUT}/${YEAST_CYCLE}_matchted_mlps.txt" > "${OUT}/sort_${YEAST_CYCLE}_matchted_mlps.txt"
awk ' $2 > 0.8 ' "${OUT}/sort_${YEAST_CYCLE}_matchted_mlps.txt" > "${OUT}/sort_${YEAST_CYCLE}_0.8_matchted_mlps.txt"
done	


