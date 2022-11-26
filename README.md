# The SAMPL9 Blind Prediction Challenges for Computational Chemistry

Challenge details, inputs, and results from the SAMPL9 series (phase) of challenges. Each individual SAMPL9 challenge may be broken up into multiple stages.

See the [SAMPL website](https://samplchallenges.github.io/) for information on the Statistical Assessement of the Modeling of Proteins and Ligands (SAMPL) series of challenges as a whole. This repository focuses specifically on the SAMPL9 series of challenges. Additionally, see the [SAMPL community](https://zenodo.org/communities/sampl/?page=1&size=20) on Zenodo for content related to the SAMPL series of challenges. If you wish to use Zenodo to post our presentation slides, datasets, and/or other SAMPL related material, and get a DOI for them, please upload to the community [here](https://zenodo.org/communities/sampl/?page=1&size=20) so that your content will be listed.

Because these files are available publicly, we have no record of who downloads them. Therefore, you should sign up for notifications. Specifically, if you want to receive updates if we uncover any problems, it is imperative that you either (a) [sign up for the SAMPL e-mail list](https://eepurl.com/dPj11j), or (b) sign up for notifications of changes to this GitHub repository (the Watch button, above); ideally you would do both.

Please note that some aspects of the SAMPL8 series of challenges are still ongoing, but as we are launching a new host-guest challenge that marks the beginning of the SAMPL9 series of challenges, so we have opened up this repository.

## Acknowledging and citing SAMPL

If you've benefitted from our work on the SAMPL series of challenges, please be sure to acknowledge our SAMPL NIH grant in any publications/presentations. This funded host-guest experiments, as well as our work organizing and administrating these challenges. You may acknowledge SAMPL by saying something like, "We appreciate the National Institutes of Health for its support of the SAMPL project via R01GM124270 to David L. Mobley (UC Irvine)."

We also ask you to cite the SAMPL dataset(s) you used. You may cite the sets by their DOI, which is available here: [![DOI](https://zenodo.org/badge/397390406.svg)](https://zenodo.org/badge/latestdoi/397390406)


Of course, we also appreciate it if you cite any overview/experimental papers relevant to the particular SAMPL challenge you participated in.

## What's here
- Host-guest challenge files for the WP6 and bCD host-guest challenges
- [Host-guest participation instructions](https://github.com/samplchallenges/SAMPL9/blob/master//host_guest_instructions.md) with information on the submissions format, etc. Submission templates are available in the subdirectories for individual host-guest systems.
- Experimental data for the WP6 and bCD challenges is available [here](https://github.com/samplchallenges/SAMPL9/blob/main/experimental_data/experimental_measurements.csv)
- Information on participation in the nanoluciferase (nanoluc) protein-ligand binding challenge for SAMPL9 and SAMPL10, in [protein_ligand/README.md](protein_ligand/README.md)

## What's coming

- Analysis of the WP6/bCD challenge
- Additional inputs for NanoLuc challenge
- Submission details and formats for NanoLuc challenge
- Info on submission of containerized methods for NanoLuc challenge

## Changes and Data Set Versions/Changelog

### Releases
- **Release 0.1** (Sept. 7, 2021, DOI [10.5281/zenodo.5485849](https://dx.doi.org/10.5281/zenodo.5485849)): WP6 challenge details, deadline, host/guest files, submission template, and submission instructions, initially added Aug. 20


### Changes not in a release
- 2021-09-29: Add link to [WP6 submission server](http://sampl-submit.us-west-1.elasticbeanstalk.com/submit/SAMPL9-HG)  
- 2021-10-30: Push deadline back to Nov. 15
- 2021-11-08: Make WP6 G4 guest optional as it contains silicon
- 2021-11-20: Add bCD challenge information, host and guest files, and submission template. Deadline Feb. 23
- 2021-11-30: Add experimental data for WP6
- 2021-12-07: Add pKa values for WP6 host in `experimental_data/WP6`
- 2021-12-08: Correct TFP `.mol2` , `.sdf`, and `.pdb` structure files in `host_guest/bCD/guest_files` directory
- 2021-12-21: Add SAMPL WP6 experimental reference at https://doi.org/10.1039/D1NJ05209H
- 2022-01-27: Add submission instructions for CD challenge
- 2022-04-06: Add CD submissions, updated usermap, and reference calculations
- 2022-04-06: Add CD raw experimental measurements provided by Gilson group
- 2022-04-06: Add generated machine readable experimental measurements for WP6 and CD
- 2022-04-25: Update analyze_hostguest.py to include CD dataset submissions.
- 2022-04-25: Edit analyze_hostguest.py updating deprecated .as_matrix() function to .to_numpy().
- 2022-04-25: Edit and update submissions.py removing stats_funcs=None parameter from plotting functions
- 2022-04-25: Add preliminary analysis of ranked methods and all methods separately for SAMPL9 WP6 and CD datasets.
- 2022-05-10: Add additional info on nanoluciferase (NanoLuc) binding prediction challenge
- 2022-06-29: Add (optional) conformers for NanoLuc challenge, as well as info on which compounds were tested as a racemic mix
- 2022-07-01: Add submission templates for NanoLuc challenge, challenge deadlines.
- 2022-07-21: Add submission link for NanoLuc Stage 1
- 2022-08-01: Correct list of compounds for which conformer generation failed, for NanoLuc Stage 1 (`tranche1_failures.txt1`)
- 2022-08-01: Add (optional) correct SMILES strings for molecules for which conformer generation had failed from both tranches. Participants may either skip these molecules or provide predictions for them.
- 2022-09-22: Extend NanoLuc Stage 1 deadline to Oct. 7
- 2022-10-10: Add tranche1 NanoLuc binders (answers to Stage 1/inputs for Stage 2)
- 2022-10-10: Add NanoLuc Stage 1 submissions
- 2022-11-16: Extend NanoLuc Stage 2 to allow for ranking rather than just IC50 prediction
- 2022-11-25: Add logP challenge initial information

## Challenge construction

### Overview

The SAMPL9 phase of challenges includes two host-guest challenges (already closed), a protein-ligand challenge on nanoluciferase, and a toluene-water logP challenge.

### The NanoLuc challenge

The SAMPL9 (and SAMPL10) protein-ligand challenges focus on binding to nanoluciferase (nanoluc) and include both a virtual screening and potency prediction component, as [further detailed here](protein_ligand/README.md). Compound identities and submission deadlines [are available](protein_ligand/README.md)], along with submission links.

## The logP challenge

A generous contribution of new experimental data makes it possible for us to run a toluene-water logP challenge, with new measurements for 16 compounds spanning a few log units in dynamic range for logP. Additional details [available under the logP directory](logP/README.md).

### The WP6 challenge

The WP6 challenge focuses on binding of WP6 to thirteen guests. Binding has been experimentally characterized, and the Isaacs group is preparing a paper for publication. Experimental results/data is now available as of 11/30/2021, [`here`](experimental_data/).

**Deadline**: The deadline for WP6 submissions was November 15, 2021. [The submission format is available here](https://github.com/samplchallenges/SAMPL9/blob/main/host_guest/WP6/WP6_submissions.txt).

**Experimental work** on [the WP6 challenge data has now been published](https://doi.org/10.1039/D1NJ05209H) in the *New Journal of Chemistry*.

### The bCD-Phenothiazine-based drugs challenge

The bCD challenge focuses on binding of five phenothiazine antipsychotic drugs to two hosts in the cyclodextrin family. Experimental binding affinity measurements have been measured, and the Gilson group is preparing a publication for this work. Experimental results/data will be available after the challenge closes.  

**Deadline**: The deadline for bCD submissions was Feb. 23, 2022. The experimental data is available here [bCD experimental data](https://github.com/samplchallenges/SAMPL9/blob/main/experimental_data/CD).

## MANIFEST
- [`host_guest/`](host_guest/): Details on host-guest challenges.
- [`experimental_data/WP6/`](experimental_data/WP6/): Experimental data for WP6 challenge provided by Isaacs on 11/30/2021. Contains: [`SAMPL9_Datasheet_20210727.cdx`](experimental_data/WP6/SAMPL9_Datasheet_20210727.cdx) - a ChemDraw file that contains smiles strings and a view of the host and guest molecules, [`host_and_guest_smiles.csv`](experimental_data/WP6/host_and_guest_smiles.csv) - a `.csv.` containing the host and guest smile strings (taken from [`SAMPL9_Datasheet_20210727.cdx`](experimental_data/WP6/SAMPL9_Datasheet_20210727.cdx)),[`SAMPL9_answersheet.docx`](experimental_data/WP6/SAMPL9_answersheet.cdx)- a MS Word document that contains the measured thermodynamic data, and [`SAMPL9_answersheet.csv`](experimental_data/WP6/SAMPL9_answersheet.csv)- a machine-readable format version of [`SAMPL9_answersheet.docx`](Sexperimental_data/WP6/SAMPL9_answersheet.cdx).
- [`experimental_measurements.X`]: Summary table of experimental data in `.csv` and `.json` formats. Includes WP6 and bCD.
- [`experimental_data/CD/`](experimental_data/CD/PhenothiazineCD-Binding-Summary-3-31-2022_BAedit.docx): Data provided by the Gilson group in `.docx` format. This document was updated/corrected on 3/31/2022 from earlier values.
- ['protein_ligand`](protein_ligand): Information on the protein-ligand challenge on nanoluciferase (nanoluc) binding, which includes library screening and prediction of potency (IC50).

## SAMPL-related
If you give a SAMPL-related talk or presentation or an analysis of its data, and are willing to share publicly, please consider posting on Zenodo and linking it to the [SAMPL Zenodo community](https://zenodo.org/communities/sampl?page=1&size=20).

Please also be sure to let David Mobley know of any SAMPL-related publications so he can include these in the [master list of all SAMPL publications](https://www.samplchallenges.org/history/allreferences/) on the [SAMPL website](samplchallenges.org)

## LICENSE

This material here is made available under a CC-BY and MIT licenses, as appropriate:

- MIT for all software/code
- CC-BY 4.0 for all other materials

In other words, we are happy to have you reuse any of the materials here for any purpose as long as proper credit/citation is given.  
