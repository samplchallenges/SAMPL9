# NCATS NanoLuc Experimental Data

This contains input information for the SAMPL9 NanoLuc Challenge, from NCATS, as detailed in the [challenge description](../README.md).

Tranche 1 of the compounds will form the basis of Phase 1 of the SAMPL9 NanoLuc challenge. A reduced selection (those which bind) will be used for Phase 2. Tranche 2 will be reserved for the SAMPL10 NanoLuc Challenge.

## Manifest
- `nanoluc_compounds_tranche1.csv`: Tranche 1 of compounds for which to predict NanoLuc binding. Contains compound IDs and isomeric SMILES. Predictions should be of whether compounds do or do not bind with detection thresholds as specified in the challenge description.
- `nanoluc_compounds_tranche2.csv`: Tranche 2 of compounds for which to predict NanoLuc binding. Contains compound IDs and isomeric SMILES. Predictions should be of whether compounds do or do not bind with detection thresholds as specified in the challenge description.
- [tranche1_conformers.sdf.gz](https://drive.google.com/file/d/1ki7Y8BgmrsbVqhsd2A6HPqnytLZDLFsk/view?usp=sharing): This file, available via Google Drive, is a 2.2 GB sdf file containing generated conformers for essentially all of the molecules in tranche1 (skipping conformer generation failures noted below)
- [tranche2_conformers.sdf.gz](https://drive.google.com/file/d/1kl6Nl9_5gq4YwKH2GQ93jXc4fkqdMtFz/view?usp=sharing): This file, available via Google Drive, is a 2.2 GB sdf file containing generated conformers for essentially all of the molecules in tranche2 (skipping conformer generation failures noted below)
- `tranche1_unspecified_stereochem.txt`: Compounds with unspecified stereochemistry in tranche 1 (for which conformer generation initially failed); these were tested as a racemic mixture, per NCATS
- `tranche2_unspecified_stereochem.txt`: Compounds with unspecified stereochemistry in tranche 2 (for which conformer generation initially failed); these were tested as a racemic mixture, per NCATS
- `tranche1_failures.txt`: Molecules in tranche1 for which we failed to generate conformers; this likely means the SMILES strings are invalid or cannot be correctly parsed.
- `tranche2_failures.txt`: Molecules in tranche2 for which we failed to generate conformers; this likely means the SMILES strings are invalid or cannot be correctly parsed.

## Additional information forthcoming
- Experimental results will be released after phases close
