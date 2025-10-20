#!/bin/bash

for i in {B,C,D,E,F,G,H,I}
do
	for j in {0,1,2}
	do
		crab resubmit -d ParkingSingleMuon${j}_2024${i}/crab_projects/crab_ParkingSingleMuon${j}_2024${i}
	done
done
