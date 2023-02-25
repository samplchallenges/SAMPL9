#load modules
import pandas as pd

#define files
ncats_exp_data = pd.read_excel('/Users/amezcum1/Desktop/secNanoLuc-Genesis-qHTS-data.xlsx')

tranche1_molecules = pd.read_csv('/Users/amezcum1/Desktop/SAMPL9/protein_ligand/NCATS_experimental_data/nanoluc_compounds_tranche1.csv')

tranche2_molecules = pd.read_csv('/Users/amezcum1/Desktop/SAMPL9/protein_ligand/NCATS_experimental_data/nanoluc_compounds_tranche2.csv')

#define each molecule ids for each tranche
tranche1_ids = []
tranche2_ids = []

for i in tranche1_molecules['Sample ID']:
    tranche1_ids.append(i)
    
for i in tranche2_molecules['Sample ID']:
    tranche2_ids.append(i)
    
#split ncats exp data
print("Splitting data...")

#tranche1
tranche1_names = []
tranche1_ic50s = []
tranche1_structures = []

#tranche2
tranche2_names = []
tranche2_ic50s = []
tranche2_structures = []

#search for tranche1 stage 2 data
for name, ic50, structure in zip(ncats_exp_data['Sample ID'], ncats_exp_data['IC50 (uM)'], ncats_exp_data['Structure']):
    if name in tranche1_ids:
        tranche1_names.append(name)
        tranche1_ic50s.append(ic50)
        tranche1_structures.append(structure)
    elif name in tranche2_ids:
        tranche2_names.append(name)
        tranche2_ic50s.append(ic50)
        tranche2_structures.append(structure)
    else:
        continue
        
#make a new df for each tranche
tranche1_dic = {'Sample ID': tranche1_names, 'IC50 (uM)': tranche1_ic50s, 'SMILES': tranche1_structures}
tranche1_ncats_exp_data = pd.DataFrame(tranche1_dic)

tranche2_dic = {'Sample ID': tranche2_names, 'IC50 (uM)': tranche2_ic50s, 'SMILES': tranche2_structures}
tranche2_ncats_exp_data = pd.DataFrame(tranche2_dic)

#drop inactives, for stage 2 results
tranche1_actives = tranche1_ncats_exp_data.dropna(subset = ['IC50 (uM)']) #contains IC50 values
tranche2_actives = tranche2_ncats_exp_data.dropna(subset = ['IC50 (uM)']) #contains IC50 values

#drop IC50 values, for stage 1 results
tranche1_binders = tranche1_actives.drop('IC50 (uM)', axis=1)
tranche2_binders = tranche2_actives.drop('IC50 (uM)', axis=1)

#save csv files
print("Saving tranche1 & tranche2 files...")
tranche1_actives.to_csv('nanoluc_binders_tranche1_with_IC50.csv', index=False)
tranche1_binders.to_csv('nanoluc_binders_tranche1.csv', index=False)

tranche2_actives.to_csv('nanoluc_binders_tranche2_with_IC50.csv', index=False)
tranche2_binders.to_csv('nanoluc_binders_tranche2.csv', index=False)

print("Done!")
