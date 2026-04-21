#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os,sys
from General import *


# In[2]:


NAME_FIELD="zincid"
SMILES_FIELD="smiles"
#binder_nonbinder_classification=lambda df: (pd.to_numeric(df["Binder"],errors="coerce").to_numpy().astype(float)>=0.5).astype(int) # 16 is the cutoff after which it qualifies as binder. Change this function to return 0s or 1s as flags for binders and non-binders
binder_nonbinder_classification=lambda df: (pd.to_numeric(df["A2aAR pKi"],errors="coerce").to_numpy().astype(float)>=5).astype(int) # 16 is the cutoff after which it qualifies as binder. Change this function to return 0s or 1s as flags for binders and non-binders
HIT_FRACTION=0.1
TOTAL_LIGS=144
SIMILARITY_CUTOFF=0.4 # How similar can hits and decoys be to each-other (in tanimoto similarity)?

os.system("rm -rf extracted_poses/*")
ligdf=pd.read_csv("data_invitro.csv")
print(ligdf.head())
lignames=list(ligdf[NAME_FIELD].to_numpy().astype(str))
ligsmi=[str(s) for s in ligdf[SMILES_FIELD].to_numpy().astype(str)]
ligclass=binder_nonbinder_classification(ligdf) #(ligdf.DispRate.to_numpy()>=16).astype(int)
for i in range(len(lignames)):
    print(lignames[i],ligclass[i])
ligdict=dict()
for ni,n in enumerate(lignames):
    ligdict[n]=ligsmi[ni]
print(len(ligdict),"ligands in total")


# In[3]:


N_HITS=-1
#N_HITS=round(TOTAL_LIGS*HIT_FRACTION)
#N_MISSES=TOTAL_LIGS-N_HITS


# In[4]:


def select(smiles,names,num,cutoff=None):
    if cutoff is None: cutoff=SIMILARITY_CUTOFF
    smiload=SequentialSMILESLoader(smiles,attach_names=names)
    simmat=getSimilarityMatrix(smiload,precompute_fps=True)
    sels=largest_dissimilar_subset(simmat,cutoff,print_freq=5000)
    _,names=smiload.drain(as_smiles=True)
    sel_n=[names[i] for i in sels]
    if len(sel_n)<num: raise ValueError(f"Not enough ligands to select with similarity cutoff of {cutoff}. Consider raising this value or picking fewer ligands")
    if num<0: return np.array(sel_n)
    else: return np.random.choice(sel_n,num,replace=False)


# In[5]:


hitlocs=np.where(ligclass)[0]
hit_smiles=[ligsmi[i] for i in hitlocs]
hit_names=[lignames[i] for i in hitlocs]
hit_sels=select(hit_smiles,hit_names,N_HITS)
N_HITS=len(hit_sels)
N_MISSES=TOTAL_LIGS-N_HITS
print("Hits:",hit_sels)
print("Total hit #:",N_HITS)
print("Targetting misses:",N_MISSES)


# In[6]:


misslocs=np.where(1-ligclass)[0]
miss_smiles=[ligsmi[i] for i in misslocs]
miss_names=[lignames[i] for i in misslocs]
miss_sels=select(miss_smiles,miss_names,N_MISSES)
print("Decoys:",miss_sels)

# In[7]:


final_keys=[]
slist=open("selected_list.smi","w")
for k in hit_sels:
    slist.write(ligdict[k]+" "+k+"_hit\n")
    final_keys.append(k)
for k in miss_sels:
    slist.write(ligdict[k]+" "+k+"_miss\n")
    final_keys.append(k)
slist.close()
print(len(final_keys),"total ligands picked")


# In[8]:


mol2prepper=SequentialMol2Loader("top500K_poses.mol2",attach_names=True)
block_data=[None]*len(final_keys)
neutrals=0
inputs=mol2prepper.getNext(1)
while inputs is not None:
    mols,names=inputs
    if mols[0].computeNetCharge()<0.2: neutrals+=1
    for ki in range(len(final_keys)):
        if block_data[ki] is not None: continue
        if final_keys[ki] in names[0]:
            block_data[ki]=mols[0]
            print(final_keys[ki],"found")
            break
    inputs=mol2prepper.getNext(1)
print("Neutrals:",neutrals)


# In[10]:

EXTRACTED_FOLDER="extracted_poses/"
hitfile=open(f"{EXTRACTED_FOLDER}/hit_data_mol2.spc","w")
unfound=0
for ki,k in enumerate(final_keys):
    blk=block_data[ki]
    if blk is None:
        print("WARN: Key",k,"not found")
        unfound+=1
        continue
    blkf=open(EXTRACTED_FOLDER+"/"+k+".mol2","w")
    blkf.write(blk.reconstructMol2Block()+"\n")
    blkf.close()
    hitfile.write(k+".mol2 "+str(ki<N_HITS)+" "+str(round(blk.computeNetCharge()))+"\n")
hitfile.close()
print("Poses not found:",unfound)
print("Poses found:",len(final_keys)-unfound)
