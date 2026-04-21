# CTMD
This github code repo provides an open-source implementation of the CTMD protocol for hit-triaging in high-throughput virtual screening campaigns. \\
If you use this method, please cite [our work](https://doi.org/10.64898/2026.02.05.703972) on which this method is based.

## Requirements
1. Linux Operating System (local or HPC are both supported)
2. [Amber MD Engine](https://ambermd.org/GetAmber.php) (pmemd24) with [PLUMED 2.9](https://www.plumed.org/download) (tested) or later.
3. Independently, AmberTools25 (for CppTraj,Antechamber,tLeap in PATH) is required for post-processing and setting up the metadynamics

**Note:** If you would use the MD scripts directly, they rely on portable Singularity Image Files (SIFs). Due to their large size, they cannot be uploaded to github, but will be provided on request

## Code organization
- All important scripts are in the base repository. The key files are:
  - `hit_data_mol2.spc`: Which contains a space-separated list of all ligands to process with CTMD. It has 3 columns (Mol2 Filename, Binder/Non-binder, Net charge of the ligand). The binder/non-binder can just be set to '1' for new screens where it is not known. This is what the file looks like
    ```
    Z2381441.mol2 1 0
    Z94372434.mol2 1 -1
    Z3596334.mol2 0 0
    Z453683233.mol2 1 +1
    ...
    ```
  - `pick_ligands.py`: Useful to pick ligands from an existing virtual screen. **You must edit this script, because it was specifically designed for retrospective studies**. Change 'binder_nonbinder_classification' on line 20 to reflect which ligands would be considered "binders" and which non-binders. The program will select a random sample of size 'TOTAL_LIGS' with a base hit-rate of 'HIT_FRACTION'. 
- **Setting up a new system:** This code is designed to work with the Mol2 file format. [OpenBabel](https://openbabel.org/docs/UseTheLibrary/PythonInstall.html) can be usedto convert between different file-formats. There should be one mol2 file for each entry in the `hit_data_mol2.spc` file with the same name. The ligands should be in the predicted binding pose (i.e. most likely docking pose) relative to the protein PDB chosen. You can then run each setup script in the following order:
  - `1_folder_structure.sh`: Sets up a folder structure with one folder for each ligand. This isolates each MD simulation from the rest and makes parallelization easier. You can either dump all the mol2 files into a folder called `extracted_poses/` and pick the molecules needed from CTMD by editing which molecules appear in the `hit_data_mol2.spc` file here. Alternatively, you can paste them all to the base folder.
  - **Fix `system_setup.leap`**: Protein and system preparation for MD can be complicated and variable. You must prepare your system before starting the procedure. `system_setup.leap` provides a template where you can combine all prepared components to build a system to simulate. All paths are relative to one subfolder of the base (i.e. a `protein.pdb` file in the base will be accessible at `../protein.pdb`).
  - `2_setup.sh`: Prepares ligand parameters with antechamber. The 'antechamber' command must be runnable and be in your 'PATH'. Make sure that the 'net_charge' entries in `hit_data_mol2.spc` are correct, or antechamber may crash. The script will inform you if this happens. Ensure that the `system_setup.leap` script doesn't crash (warnings are fine, but the last line must say '0 errors').
  - `./start_MD.sh \[foldername\]`: This will start the MD for each folder. Note that if you are using this script directly, you will need our singularity image file (see above). We are currently working on making it accessible, but due to its large size, it cannot be uploaded to GitHub. Alternatively, you can edit the the `amber_MD_template/submit_MD.sh` script to instead directly run locally if you have Pmemd24 installed locally.
