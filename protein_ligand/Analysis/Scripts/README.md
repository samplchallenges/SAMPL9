# SAMPL9 NanoLuc analysis scripts

## Manifest
- `analyze_stage1.py`: Master analysis script, for stage1 of the SAMPL9 NanoLuc Challenge.
- `get_usermap.py`: Should be run before `analyze_stage1.py`; obtains and stores a submission map which lists submission IDs and corresponding file names.
- `pkganalysis`: Utility classes/functions to parse SAMPL submissions and for statistical analysis.
- `hits_verification.csv`: file generated from nanoluc_binders_tranch1.csv, to include all molecules 
