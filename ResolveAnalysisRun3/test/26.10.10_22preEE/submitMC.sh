#!/bin/bash

name_list=()
data_list=()

while read -r name data
do
        name_list+=("$name")
        data_list+=("$data")
done < inputdataMC.dat

for i in {0..35}
do
        echo "${name_list[$i]}"
        mkdir -p "crab/${name_list[$i]}"

        if (( i >= 21 && i <= 35 )); then
                submit_file="submitSig_crab.py"
        else
                submit_file="submitMC_crab.py"
        fi

        sed -e "s/test_1/${name_list[$i]}/g" \
            "$submit_file" > "./crab/${name_list[$i]}/submit_crab_tmp.py"

        sed -e "s?dataset?${data_list[$i]}?g" \
            "./crab/${name_list[$i]}/submit_crab_tmp.py" \
            > "./crab/${name_list[$i]}/${submit_file}"

        cp runResolvedNtuplizer_run3MC_cfg.py "./crab/${name_list[$i]}"

        (
                cd "./crab/${name_list[$i]}" || exit 1
                crab submit "$submit_file"
        )
done
