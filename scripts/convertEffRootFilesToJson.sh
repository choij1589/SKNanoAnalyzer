#!/bin/bash
ERA=$1

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/MUO/root/efficiency_TopHNT.root \
  --inputHist data:sim:sf --outputName eff_data:eff_mc:sf \
  --csetName "Measured efficiency for TopHNT ID" \
  --outputPath ${SKNANO_DATA}/${ERA}/MUO/efficiency_TopHNT.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/MUO/root/efficiency_Mu17Leg1.root \
  --inputHist data:sim:sf --outputName eff_data:eff_mc:sf \
  --csetName "Measured efficiency for Mu17Leg1 Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/MUO/efficiency_Mu17Leg1.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/MUO/root/efficiency_Mu8Leg2.root \
  --inputHist data:sim:sf --outputName eff_data:eff_mc:sf \
  --csetName "Measured efficiency for Mu8Leg2 Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/MUO/efficiency_Mu8Leg2.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/MUO/root/efficiency_Mu8El23_Mu8Leg.root \
  --inputHist Mu8El23_Data:Mu8El23_MC --outputName eff_data:eff_mc \
  --csetName "Measured efficiency for Mu8El23_Mu8Leg Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/MUO/efficiency_Mu8El23_Mu8Leg.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/MUO/root/efficiency_Mu23El12_Mu23Leg.root \
  --inputHist Mu23El12_Data:Mu23El12_MC --outputName eff_data:eff_mc \
  --csetName "Measured efficiency for Mu23El12_Mu23Leg Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/MUO/efficiency_Mu23El12_Mu23Leg.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/EGM/root/efficiency_TopHNT.root \
  --inputHist data:sim:sf --outputName eff_data:eff_mc:sf \
  --csetName "Measured efficiency for TopHNT ID" \
  --xaxis sceta --yaxis pt \
  --outputPath ${SKNANO_DATA}/${ERA}/EGM/efficiency_TopHNT.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/EGM/root/efficiency_Mu8El23_El23Leg.root \
  --inputHist Mu8El23_Data:Mu8El23_MC --outputName eff_data:eff_mc \
  --xaxis abssceta --yaxis pt \
  --csetName "Measured efficiency for Mu8El23_El23Leg Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/EGM/efficiency_Mu8El23_El23Leg.json

./scripts/convertRootToJson.py \
  --inputFile ${SKNANO_DATA}/${ERA}/EGM/root/efficiency_Mu23El12_El12Leg.root \
  --inputHist Mu23El12_Data:Mu23El12_MC --outputName eff_data:eff_mc \
  --xaxis abssceta --yaxis pt \
  --csetName "Measured efficiency for Mu23El12_El12Leg Trigger path" \
  --outputPath ${SKNANO_DATA}/${ERA}/EGM/efficiency_Mu23El12_El12Leg.json

