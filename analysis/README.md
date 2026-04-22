# Analysis scripts
Attached are the analysis scripts used for our work. These scripts use a bootstrapping method to subselect ligands (where possible) to generate enrichment plots with uncertainties based on which binders were picked.

## Requirements
In order to run these, you will need to have the following packages installed
- numpy
- matplotlib
- pandas
- rdkit
- tqdm
Most of these are available via pip (i.e. pip install numpy matplotlib pandas rdkit tqdm)

## Using the scripts
These scripts are designed to use our specific file-structure and work for each of our 4 systems. **Only for JAK2**, the source of the data was different, and ABFEP data was also available.\\
To get our analysis results for this system, use: `analysis_withABFEP.py` \\
For all other systems, use: `analysis_noABFEP.py [AmpC/Alpha2AR/CB1R]`

## Generic new problem
If you are working on a new system with no known binder information (i.e. for a prospective run), use `./analysis_unknown.py`
