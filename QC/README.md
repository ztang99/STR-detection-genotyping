# Quality Control

Scripts used for conducting sequencing quality controls.

- `IPNSamtoolsCheck.sh`: Quality check of BAM files using samtools
- `ExtractInfo4Table.sh`: Extract QC metrics from samtools output
- `ethnicity_pred_gnomad_cont.py`: Predict sample ancestry using gnomAD continental probabilities
- `snvstory_ethnicity_check.sh`: Run SNVstory for ancestry determination
- `snvstory_selectVariant.sh`: Select variants from gVCF files for SNVstory