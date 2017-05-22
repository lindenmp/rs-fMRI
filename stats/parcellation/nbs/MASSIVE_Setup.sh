#!/bin/env bash

for WhichProject in OCDPG UCLA; do
	for WhichSplit in Diagnostic Motion; do
		for WhichParc in Gordon Power; do
			for WhichNoise in 6P 6P+2P 6P+2P+GSR 24P 24P+8P 24P+8P+4GSR 24P+aCC 24P+aCC+4GSR 24P+aCC50 24P+aCC50+4GSR 12P+aCC 12P+aCC50 sICA-AROMA+2P sICA-AROMA+2P+GSR sICA-AROMA+8P sICA-AROMA+8P+4GSR; do
				echo ${WhichProject} ${WhichSplit} ${WhichParc} ${WhichNoise}
				sbatch MASSIVE_Job.sh ${WhichProject} ${WhichSplit} ${WhichParc} ${WhichNoise} 0.05 10k_NBS_tDOF false
				sleep 1
			done
		done
	done
done

echo -e "\n***FINISHED***\n"
