#!/bin/bash

##############################################################################

# Master script to run the entire pipeline.

## author: Zitian Tang
## contact: tang.zitian@wustl.edu

##############################################################################

SUBNAME="$1"
PROJECT_DIR="."

cd $PROJECT_DIR || exit 1

## 1_EHdn_GenerateStrProfile
bash STR_detection_pipeline/1_EHdn_GenerateStrProfile.sh input/cases.txt cases output
bash STR_detection_pipeline/1_EHdn_GenerateStrProfile.sh input/controls.txt controls output


## 2_EHdn_GenerateManifestFile
bash STR_detection_pipeline/2_EHdn_GenerateManifestFile.sh $PROJECT_DIR $SUBNAME \
    input/cases.txt \
    output/EHdn/$SUBNAME/EHdn_cases_str-profiles \
    input/controls.txt \
    output/EHdn/$SUBNAME/EHdn_controls_str-profiles

## 3_EHdn_RunAnnotEHdn
bash STR_detection_pipeline/3_EHdn_RunAnnotEHdn.sh $PROJECT_DIR $SUBNAME

## 4_EH_RunEH
bash STR_detection_pipeline/4_EH_RunEH.sh $PROJECT_DIR $SUBNAME

## 4.5_CombineEHResult.sh
bash STR_detection_pipeline/4.5_CombineEHResult.sh $PROJECT_DIR $SUBNAME

## 5_RunBLAT
bash STR_detection_pipeline/5_RunBLAT.sh $PROJECT_DIR $SUBNAME input/cases.txt input/controls.txt input/roi.bed

## 6_BuildDatabase
bash STR_detection_pipeline/6_BuildDatabase.sh $PROJECT_DIR $SUBNAME input/roi.bed

## 7_QueryDatabase
bash STR_detection_pipeline/7_QueryDatabase.sh $PROJECT_DIR $SUBNAME input/roi.bed