#!/bin/bash

name_list=()
data_list=()

while read name data
do
	name_list+=("$name")
	data_list+=("$data")
done < inputdata2.dat

for i in {0..23}
do
	echo ${name_list[$i]}
	#mkdir -p crab/${name_list[$i]}
	#sed -e "s/test_1/${name_list[$i]}/g" submit2_crab.py > ./crab/${name_list[$i]}/submit_crab_tmp.py
	#sed -e "s?dataset?${data_list[$i]}?g" ./crab/${name_list[$i]}/submit_crab_tmp.py > ./crab/${name_list[$i]}/submit_crab.py
	#cp runMergedLeptonIDImageMC_run3_cfg.py ./crab/${name_list[$i]}
	cd ./crab/${name_list[$i]}
	#crab submit submit_crab.py
	#crab status -d crab_projects/crab_Image_match_${name_list[$i]} 
	#crab kill -d crab_projects/crab_Image_${name_list[$i]} 
	crab resubmit -d crab_projects/crab_Image_match_${name_list[$i]} 
	cd ../..
done
