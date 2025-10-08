#!/bin/bash

ERAs=("2016preVFP" "2016postVFP" "2017" "2018" "2022" "2022EE" "2023" "2023BPix")
OBJs=("muon" "electron")

for era in "${ERAs[@]}"; do
    for obj in "${OBJs[@]}"; do
        echo "Processing era: $era, object: $obj"
        ./scripts/submission/MeasFakeRate.sh $era $obj
    done
done
