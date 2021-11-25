# Guests for the SAMPL9 bCD Challenge

Provided are the guests for this challenge in PDB, Tripos MOL2 and SDF file formats. The guests are codenamed from `TDZ`, `TFP`, `PMZ`, `PMT`, and `CPZ`..

## What's here

- `guest_smiles.txt`: Source file hand-generated from Gilson-provided `SAMPL9-Host.cdx` by copying-and-pasting SMILES (by Martin Amezcua, November 19, 2021).
- `input_maker.ipynb`: The jupyter notebook used to generate the PDB, MOL2 and SDF files for each guest using OpenEye toolkits and the SMILES strings and codenames found in `guest_smiles.txt`. For some guests, random stereoisomers were chosen; steroechemistry of the nitrogen center may require care.
- `.mol2`, `.sdf` and `.pdb` files for all guests: These were auto-generated from the SMILES via `input_maker.ipynb` using the OpenEye toolkits.
