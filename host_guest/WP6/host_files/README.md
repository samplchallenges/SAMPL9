# WP6 host files

The WP6 files here were generated from the `WP6_host_smiles.txt` SMILES file created manually. 

## Manifest
- `WP6_host_smiles.txt`: SMILES of the WP6 host. (Manually created by Martin Amezcua on August 18, 2021)
- `WP6.sdf`: created from SMILES WP6_host_smiles using openeye toolkit
- `WP6.mol2`: mol2 format WP6 file generated from WP6_host_smiles.txt. AM1BCC charges were assigned using `Assign_Charge.ipynb`.
- `WP6.pdb`: PDB format WP6 file
- `Assign_Charge.ipynb`: Jupyter notebook used to assign AM1BCC charges (via Openeye) to `WP6.mol2` host file.  
