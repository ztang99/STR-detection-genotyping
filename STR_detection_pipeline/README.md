# STR Detection Pipeline

Scripts for the actual detection pipeline. All helper scripts called by bash scripts are stored inside `python_scripts/`.

`00_RunAll.sh`: Master script to run all steps at once (for local environment with all necessary software installed). Individual steps shown below can be run separately as needed.

`1_EHdn_GenerateStrProfile.sh`: Generate STR profiles using ExpansionHunterDenovo

`2_EHdn_GenerateManifestFile.sh`: Generate manifest file for ExpansionHunterDenovo

`3_EHdn_RunAnnotEHdn.sh`: Run and annotate ExpansionHunterDenovo results\
Calls helper script: `python_scripts/wdl_filter_ehdn_results.py`

`4_EH_RunEH.sh`: Run ExpansionHunter on detected STR regions\
Calls helper script: `python_scripts/wdl_IPN_generate_EHcatalog.py`

`5_CombineEHResult.sh`: Combine ExpansionHunter results\
Calls helper script: `python_scripts/wdl_filter_eh_vcfs.py`

`6_CombineEhdnEH.sh`: Obtain consensus calls from ExpansionHunterDenovo and ExpansionHunter results\
Calls helper script: `python_scripts/wdl_combine_ehdn_eh.py`

`7_RunBLAT.sh`: Run BLAT alignment of STR regions\
Calls helper script: `python_scripts/wdl_query_STR_db.py`

`8_BuildDatabase.sh`: Build STR sequence database\
Calls helper script: `python_scripts/wdl_addBlatResult2db.py`

`9_QueryDatabase.sh`: Query STR database\
Calls helper script: `python_scripts/wdl_query_STR_db.py`

`10_UnsupervisedIPNFinalStep.R`: R script for unsupervised clustering of repeat expansions
