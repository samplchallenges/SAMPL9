# NCATS NanoLuc Experimental Data

This contains input information for the SAMPL9 NanoLuc Challenge, from NCATS, as detailed in the [challenge description](../README.md).

Tranche 1 of the compounds will form the basis of Phase 1 of the SAMPL9 NanoLuc challenge. A reduced selection (those which bind) will be used for Phase 2. Tranche 2 will be reserved for the SAMPL10 NanoLuc Challenge.

Some SMILES strings provided by NCATS were apparently incorrect and could not be parsed; compounds with these SMILES strings are reflected in `trancheN_failures.txt`. Later (May 2022), NCATS provided sdf files containing conformers for these molecules, and then we generated updated SMILES strings. Because of the lateness of this data, we have not updated the `nanoluc_compounds_trancheN.csv` files with these SMILES strings, but instead provided them in `nanoluc_compounds_trancheN_failure_updates.csv`. Participants can choose to include these compounds or skip them.

For Stage 2, the `nanoluc_binders_tranche1.csv` file contains compound IDs for binders, as well as original NCATS-provided SMILES strings for these. SDFs and updated/corrected SMILES strings for failures may be obtained from other files in this repository.

## Manifest
- `nanoluc_compounds_tranche1.csv`: Tranche 1 of compounds for which to predict NanoLuc binding. Contains compound IDs and isomeric SMILES. Predictions should be of whether compounds do or do not bind with detection thresholds as specified in the challenge description.
- `nanoluc_compounds_tranche2.csv`: Tranche 2 of compounds for which to predict NanoLuc binding. Contains compound IDs and isomeric SMILES. Predictions should be of whether compounds do or do not bind with detection thresholds as specified in the challenge description.
- [tranche1_conformers.sdf.gz](https://drive.google.com/file/d/1ki7Y8BgmrsbVqhsd2A6HPqnytLZDLFsk/view?usp=sharing): This file, available via Google Drive, is a 2.2 GB sdf file containing generated conformers for essentially all of the molecules in tranche1 (skipping conformer generation failures noted below)
- [tranche2_conformers.sdf.gz](https://drive.google.com/file/d/1kl6Nl9_5gq4YwKH2GQ93jXc4fkqdMtFz/view?usp=sharing): This file, available via Google Drive, is a 2.2 GB sdf file containing generated conformers for essentially all of the molecules in tranche2 (skipping conformer generation failures noted below)
- `tranche1_unspecified_stereochem.txt`: Compounds with unspecified stereochemistry in tranche 1 (for which conformer generation initially failed); these were tested as a racemic mixture, per NCATS
- `tranche2_unspecified_stereochem.txt`: Compounds with unspecified stereochemistry in tranche 2 (for which conformer generation initially failed); these were tested as a racemic mixture, per NCATS
- `tranche1_failures.txt`: Molecules in tranche1 for which we failed to generate conformers; this likely means the SMILES strings are invalid or cannot be correctly parsed. Corrected 2022-08-01. Alternate input for these failures provided in a separate file.
- `tranche2_failures.txt`: Molecules in tranche2 for which we failed to generate conformers; this likely means the SMILES strings are invalid or cannot be correctly parsed. Alternate input for these failures provided in a separate file.
- `nanoluc_compounds_tranche1_failure_updates.csv` and `nanoluc_compounds_tranche2_failure_updates.csv`: Updated (corrected) SMILES strings for molecules for which conformer generation failed
- `nanoluc_binders_tranche1.csv`: NanoLuc binders from tranche1 (compound ID and SMILES); provides answers for Stage 1, and inputs for Stage 2. SMILES strings are original NCATS SMILES strings, and corrected SMILES from the above files can be used if failures are encountered. 
- `nanoluc_binders_tranche1_with_IC50.csv`: NanoLuc binders IC50 measurements from tranche1 (compound ID, IC50, and SMILES). SMILES strings are original NCATS SMILES strings, and corrected SMILES from the above files can be used if failures are encountered.
- `Scripts`: Contains scripts used to process and/or split NanoLuc experimental data. 

## Additional information forthcoming
- Experimental results will be released after phases close
