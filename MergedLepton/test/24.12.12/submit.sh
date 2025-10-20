#!/bin/bash

prompt=''
list=()

for i in {B,C,D,E,F,G,H,I}
do
	if [ "$i" == 'E' -o "$i" == 'I' ]; then
		prompt='v2'
	else 
		prompt='v1'
	fi	
	if [ "$i" == 'G' -o "$i" == 'H' -o "$i" == 'I' ]; then
		list=`seq 0 11`
	else 
		list=`seq 0 7`
	fi	
	
	echo $i_$prompt
	for j in $list
	do
		echo $j
		#mkdir ParkingSingleMuon${j}_2024${i}
		#sed -e "s/test_1/ParkingSingleMuon${j}_2024${i}/g" submit_crab.py > ParkingSingleMuon${j}_2024${i}/submit_crab_tmp.py
		#sed -e "s?dataset?/ParkingSingleMuon${j}/Run2024${i}-PromptReco-$prompt/MINIAOD?g" ParkingSingleMuon${j}_2024${i}/submit_crab_tmp.py > ParkingSingleMuon${j}_2024${i}/submit_crab.py
		#cp runMergedLeptonIDJpsiAnalyzerData_run3_2024_cfg.py ParkingSingleMuon${j}_2024${i}/
		cd ParkingSingleMuon${j}_2024${i}/
		crab status -d crab_projects/crab_ParkingSingleMuon${j}_2024${i}
		#source /cvmfs/cms.cern.ch/crab3/crab.sh pre
		#crab submit submit_crab.py
		cd ..
	done
done
