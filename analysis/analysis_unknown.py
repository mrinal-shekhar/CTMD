#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import os,sys
import pandas as pd
from rdkit.ML.Scoring import Scoring
from default_params import *
import tqdm
from sklearn.linear_model import LinearRegression
plt.rcParams.update({"text.usetex": True, "figure.dpi": 300, "figure.constrained_layout.use": True})

#get_ipython().run_line_magic('run', '../DOCK/DOCK_OutputParser.ipynb')


# In[2]:


curfol=sys.argv[1] # Pick folder name

MIN_COLVARS=5 # Minimum number of COLVAR files (replicas) to be present for the ligand to be considered
BENCHMARK_LIGAND=None #"L95" # If any known ligand is part of the analysis, it will benchmark other scores relative to this (otherwise, choose None)
SCORE_TOLERENCE=2.49 # 1 kT
LIGAND_SMILES="/path/to/smiles.smi" # List of all ligands (must include all used in CTMD, can include many more)
LIGAND_DATA_FILE="/ath/to/saved/stats.npz" # Stats will be saved here for repeated runs. For the first time, they are automatically computed
LOW_RAM=True

hit_data_spc=np.loadtxt(f"{curfol}/hit_data_mol2.spc",dtype=str)
allowed_names=[v.replace(".mol2","") for v in hit_data_spc[:,0]]
ligand_charges=hit_data_spc[:,-1].astype(float)
print(len(allowed_names),"ligands in Hit Data")

MD_folds=sorted([fol for fol in os.listdir(curfol) if os.path.isdir(curfol+"/"+fol+"/") and fol.replace("_out","") in allowed_names])
assert len(MD_folds)==len(allowed_names), "Not all folders were found?"


# In[3]:


def processCOLVAR(cvfile,cutoff=0.35,commit_frames=2,frame_limit=-1):
    cvdata=np.loadtxt(cvfile,usecols=(0,1,-1))[:frame_limit] # Time, CV, c(t)
    if np.prod(cvdata.shape)<3: return np.nan,0
    if len(cvdata)<commit_frames:
        try: return cvdata[-1,-1],cvdata[-1,0]-cvdata[0,0]
        except: return np.nan,0

    detach=(cvdata[:,1]>cutoff)
    convfilt=np.ones(commit_frames,dtype=float)
    detach=np.convolve(detach,convfilt,mode="valid")
    retidx=np.argmax(detach)
    detqual=(detach[retidx]/commit_frames)
    if detqual<0.99: # Trajectory went to end
        retidx=len(cvdata)-1
    return cvdata[retidx,-1],cvdata[retidx,0]-cvdata[0,0] #cvdata,retidx,cvdata[retidx,-1]
def process_folder(foldpath,bootstrap_iters=250,n_repeats=3):
    colvars=sorted([f for f in os.listdir(foldpath) if f.startswith("COLVAR_")])
    #print(len(colvars),"COLVAR files loaded")
    cts=[]
    rts=[]
    for i in range(len(colvars)):
        try:
            ct,rt=processCOLVAR(foldpath+"/"+colvars[i],cutoff=CUTOFF,commit_frames=COMMIT_FRAMES,frame_limit=FRAME_LIMIT)
        except ValueError:
            print("Error processing folder:",foldpath,colvars[i])
            raise ValueError()
        cts.append(ct)
        rts.append(rt)
    cts,rts=np.array(cts,dtype=float),np.array(rts,dtype=float)
    if len(cts)<n_repeats or len(cts)<MIN_COLVARS: return np.nan,np.nan,0,0
    ct_mins=[]
    rt_mins=[]
    for _ in range(bootstrap_iters):
        idxsel=np.random.choice(np.arange(len(cts)),n_repeats,replace=False)
        ctsel=cts[idxsel]
        rtsel=rts[idxsel]
        ctidx=np.argmin(ctsel)
        ct_mins.append(ctsel[ctidx])
        rt_mins.append(rtsel[ctidx])
    return np.mean(ct_mins),np.mean(rt_mins),np.std(ct_mins),np.std(rt_mins)


# In[4]:


ligand_scores=[]
for idx in tqdm.tqdm(range(len(MD_folds)),ncols=100):
    ligand_scores.append(process_folder(curfol+"/"+MD_folds[idx])[0])
ligand_scores=np.array(ligand_scores)
ligand_names=np.array(MD_folds,dtype=str)

srt=np.argsort(-ligand_scores)
ligand_names=ligand_names[srt]
ligand_scores=ligand_scores[srt]

summary_dict=dict()
for ni,n in enumerate(ligand_names):
    summary_dict[n]=ligand_scores[ni]

if BENCHMARK_LIGAND is not None:
    if np.isnan(summary_dict[BENCHMARK_LIGAND]):
        print(f"As {BENCHMARK_LIGAND} has no score, using a default value")
        summary_dict[BENCHMARK_LIGAND]=109.69


# ### Writing Ligand Data
# Properties
if os.path.exists(LIGAND_DATA_FILE):
    raw_data=np.load(LIGAND_DATA_FILE,allow_pickle=True)
    smiles_string=raw_data["smiles_string"].item()
    if not LOW_RAM:
        rot_bonds=raw_data["rot_bonds"].item()
        mol_wt=raw_data["mol_wt"].item()
        cons_rot_bonds=raw_data["cons_rot_bonds"].item()
    else: del raw_data
else:
    # For now, just the SMILES will do
    smiles_string=dict()
    smif=open(LIGAND_SMILES,"r")
    for l in smif:
        l=l.strip().split()
        smiles_string[l[1].strip()]=l[0].strip()


# In[7]:
os.system(f'rm -rf {curfol}/results_test')
os.system(f'mkdir -p {curfol}/results_test')

plt.title(curfol.replace("/","").replace("_"," "),fontsize=28)
plt.hist(ligand_scores,bins=24)
if BENCHMARK_LIGAND is not None:
    plt.axvline(x=summary_dict[BENCHMARK_LIGAND],c="red",linestyle="--")
plt.xlabel("CTMD Score (higher is better)",fontsize=21)
plt.ylabel("No. of ligands",fontsize=21)
plt.savefig(f"{curfol}/results_test/ScoreDistribution.png")


# In[10]:


if BENCHMARK_LIGAND is not None:
    top_ligands=[k for k in summary_dict if summary_dict[k]>summary_dict[BENCHMARK_LIGAND]-SCORE_TOLERENCE and k!=BENCHMARK_LIGAND]
    print(len(top_ligands),"ligands score above "+BENCHMARK_LIGAND)
else:
    top_ligands=list(summary_dict.keys())[:int(len(summary_dict)*0.3)] # Take top 30% by default
    print(len(top_ligands),"ligands chosen")

lfo=open(curfol+"/results_test/selected_ligands.smi","w")
for n in top_ligands:
    if n in smiles_string:
        lfo.write(smiles_string[n]+" "+n+"\n")
lfo.close()

for ln in top_ligands:
    os.system(f"cp {curfol}/{ln}.mol2 {curfol}/results_test/")
np.save(f"{curfol}/results_test/CTMD_scores.npy",summary_dict)


# In[12]:

# Draw the ligands (useful in jupyter notebooks)
'''
from rdkit import Chem
from rdkit.Chem import Draw
Draw.MolDrawOptions.legendFontSize.setter(21)
ligand_objs=[Chem.MolFromSmiles(smiles_string[n]) for n in top_ligands if n in smiles_string]

#plt.figure(figsize=(16,12))
img=Chem.Draw.MolsToGridImage(
    ligand_objs[:],
    legends=top_ligands[:],subImgSize=(350,200),molsPerRow=3
)
img
'''


# In[23]:


if BENCHMARK_LIGAND is not None:
    N_BARS=30
else:
    N_BARS=len(top_ligands)
if N_BARS>len(summary_dict): N_BARS=len(summary_dict)

from matplotlib.font_manager import FontProperties
plt.rcParams.update({"font.family": "Times New Roman"})
#fontp = FontProperties(family='Century', weight='bold', size=16, math_fontfamily='dejavusans')

fig,ax=plt.subplots()
allbars=ax.barh(np.arange(N_BARS)[::-1],list(summary_dict.values())[:N_BARS])
#allbars[0].set_color("red")
barcol="green"
lig_FE=0.0
for i in range(N_BARS):
    if ligand_names[i]!=BENCHMARK_LIGAND:
        allbars[i].set_color(barcol)
        continue
    allbars[i].set_color("yellow")
    barcol="orange"
ax.set_yticks(np.arange(N_BARS)[::-1])
ax.set_yticklabels(list(summary_dict.keys())[:N_BARS])
ax.set_xlabel("$C(t)$ Score (higher is better)", fontsize=16)
plt.savefig(f"{curfol}/results_test/CTMD_bars.svg")
plt.show()

print("Selected top ligands by CTMD:")
for k in list(summary_dict.keys())[:N_BARS]:
    print("\t",k)
print("",flush=True)




