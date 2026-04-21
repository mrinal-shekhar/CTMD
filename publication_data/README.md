# Data used in our Publication
Please refer to [our work](https://doi.org/10.64898/2026.02.05.703972) for details on the protocol. These folders contain an organized copy of our simulation results.

## Folder organization
Each of the 4 systems we tested is in its own folder. This gives us 4 folders. Each folder here has all the common files:
- `hit_data_mol2.spc`
- Relevant PDBs
- Ligand poses (as mol2 files, 1 per ligand)
- One folder per ligand, containing the COLVAR files


## Ligand Folders and COLVAR files
Each ligand is simulated in one folder bearing its name. 10 replicas of metadynamics were performed for each ligand.
The COLVAR files are the main results from these metadynamics simulations. The key columns are:
1. time
2. d (RMSD of ligand to 'ref.pdb')
3. bias (Total bias at that bin - not used)
4. Cumulative bias C(t)

Each ligand folder should have around 10 copies (10 replica simulations)

## PDB files
We provide the following PDB files for each system:
- Prepared amber system as 'system.pdb'
- Just the protein (used in tleap) as 'protein.pdb'
- Reference after equilbriation (used as reference for ligand RMSD calculation) as 'ref.pdb'
- The plumed file used for the metadynamics simulation is common to all, and the template can be accessed at `amber_MD_template/plumed_template.dat` in the repo root

## hit_data_mol2.spc
This is the general data file that contains a space-separated list of all ligands to process with CTMD. It has 3 columns (Mol2 Filename, Binder/Non-binder, Net charge of the ligand). The binder/non-binder can just be set to '1' for new screens where it is not known. This is what the file looks like
```
Z2381441.mol2 1 0
Z94372434.mol2 1 -1
Z3596334.mol2 0 0
Z453683233.mol2 1 +1
...
```
