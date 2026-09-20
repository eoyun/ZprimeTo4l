#!/bin/bash

name_list=()
data_list=()

while read name data
do
	name_list+=("$name")
	data_list+=("$data")
done < inputdataMC.dat

for i in {0..1}
do
	echo ${name_list[$i]}
	mkdir -p crab/${name_list[$i]}
	sed -e "s/test_1/${name_list[$i]}/g" submitMC_crab.py > ./crab/${name_list[$i]}/submit_crab_tmp.py
	sed -e "s?dataset?${data_list[$i]}?g" ./crab/${name_list[$i]}/submit_crab_tmp.py > ./crab/${name_list[$i]}/submitMC_crab.py
	cp runMuonJpsi_run3MC_cfg.py ./crab/${name_list[$i]}
	cd ./crab/${name_list[$i]}
	crab submit submitMC_crab.py
	#crab status -d crab_projects/crab_mergedMuon_${name_list[$i]} 
	#crab resubmit -d crab_projects/crab_mergedMuon_${name_list[$i]} 
	cd ../..
done
