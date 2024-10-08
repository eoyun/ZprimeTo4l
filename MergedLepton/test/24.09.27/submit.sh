#!/bin/bash

for i in {C,D,E}
do
	for j in {0,1,2}
	do
		cd ParkingSingleMuon${j}_2022${i}/
		crab submit submit_crab.py
		cd ..
	done
done
