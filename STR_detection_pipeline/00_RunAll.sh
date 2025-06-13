#!/bin/bash

PROJECT_DIR="RunIPNonly_20250204"
SUBNAME="All"

cd $PROJECT_DIR || exit 1

## 1_EHdn_GenerateStrProfile
bash AllScripts/1_EHdn_GenerateStrProfile.sh input/cases.txt cases output
bash AllScripts/1_EHdn_GenerateStrProfile.sh input/controls.txt controls output


## 2_EHdn_GenerateManifestFile
bash AllScripts/2_EHdn_GenerateManifestFile.sh $PROJECT_DIR $SUBNAME input/cases.txt output/EHdn/EHdn_cases_str-profiles input/controls.txt output/EHdn/EHdn_controls_str-profiles

## 3_EHdn_RunAnnotEHdn
bash AllScripts/3_EHdn_RunAnnotEHdn.sh $PROJECT_DIR $SUBNAME

## 4_EH_RunEH
bash AllScripts/4_EH_RunEH.sh $PROJECT_DIR $SUBNAME

## 4.5_CombineEHResult.sh
bash AllScripts/4.5_CombineEHResult.sh $PROJECT_DIR $SUBNAME

## 5_RunBLAT
bash AllScripts/5_RunBLAT.sh $PROJECT_DIR $SUBNAME input/cases.txt input/controls.txt input/rfc1_roi.bed

## 6_BuildDatabase
bash AllScripts/6_BuildDatabase.sh $PROJECT_DIR $SUBNAME input/rfc1_roi.bed

## 7_QueryDatabase
bash AllScripts/7_QueryDatabase.sh $PROJECT_DIR $SUBNAME input/rfc1_roi.bed