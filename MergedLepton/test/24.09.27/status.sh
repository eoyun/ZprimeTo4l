#!/bin/bash

for i in {C,D,E}
do
	for j in {0,1,2}
	do
		crab status -d ParkingSingleMuon${j}_2022${i}/crab_projects/crab_240927_ParkingSingleMuon${j}_2022${i}
	done
done
