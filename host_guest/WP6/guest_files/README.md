# Guests for the SAMPL8 WP6 Challenge

Provided are the guests for this challenge in PDB, Tripos MOL2 and SDF file formats. The guests are codenamed from `G1` through `G13`.

## What's here

- `WP6_guest_smiles.txt`: Source file hand-generated from Isaacs-provided `SAMPL9 Datasheet 20210727.cdx` by copying-and-pasting SMILES (by Martin Amezcua, August 18, 2021).
- `input_maker.ipynb`: The jupyter notebook used to generate the PDB, MOL2 and SDF files for each guest using OpenEye toolkits and the SMILES strings and codenames found in `WP6_guest_smiles.txt`. For some guests, random stereoisomers were chosen; steroechemistry of the nitrogen center may require care.
- `.mol2`, `.sdf` and `.pdb` files for all guests: These were auto-generated from the SMILES via `input_maker.ipynb` using the OpenEye toolkits.
