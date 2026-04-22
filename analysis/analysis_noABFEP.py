#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import pandas as pd
from common_analysis import *
import tqdm
from rdkit.ML.Scoring import Scoring
from default_params import *

TOPK_MIN=3
N_BOOTSTRAP=360
DOCK_LABEL="DOCK 3.7" # Glide XP for JAK2 only
SYSTEM_NAME="AmpC"
if len(sys.argv)>1: SYSTEM_NAME=sys.argv[1].strip()

SETUP_FILE = f"publication_data/{SYSTEM_NAME}/hit_data_mol2.spc" # Defines path to base directory and the list of molecules
BASE_DIR=SETUP_FILE[:SETUP_FILE.rfind("/")+1]
DOCK_SCORES=BASE_DIR+"/docking_scores.spc"
BOLTZ2_FILE=f"publication_data/{SYSTEM_NAME}/{SYSTEM_NAME}_boltz_prediction.csv"
# Mutation tests only available for JAK2 and Alpha2AR
#BOLTZ2_FILE=f"publication_data/{SYSTEM_NAME}/{SYSTEM_NAME}_phe_mut_boltz_prediction.csv" # Active site mutated to PHE
#BOLTZ2_FILE=f"publication_data/{SYSTEM_NAME}/{SYSTEM_NAME}_ala_mut_boltz_prediction.csv" # Active site mutated to ALA
#BOLTZ2_FILE=f"publication_data/{SYSTEM_NAME}/{SYSTEM_NAME}_flip_mut_boltz_prediction.csv" # For boltz2 with mutations (charged residues flipped charges)

def generate_dict(names,scores,setupfile_test=None):
    ret_dict=dict();
    for ni,n in enumerate(names): ret_dict[n]=scores[ni]
    if setupfile_test is not None:
        print("*** Missing ligand list (should be empty):")
        for k in ret_dict:
            if k not in setupfile_test.names:
                print("\t",k)
                raise ValueError("List not empty! Some ligand in scores list is not part of the original")
        print("*** List ended")
    return ret_dict
setupfile=CTMDDatabase(SETUP_FILE)
print(np.sum(setupfile.binder_flag),"binders in dataset of a total of",len(setupfile.names),"ligands",flush=True)
N_HIT=min(np.sum(setupfile.binder_flag),5) # At-most take 5 binders. If fewer are available, use all for analysis

# Compute C(t) numbers
selected_ligands=np.concatenate((setupfile.sampleByClass(N_HIT,binder=True),setupfile.sampleByClass(binder=False)))
Ntot=len(selected_ligands)
BASE_HIT_RATE=(N_HIT/Ntot)
hit_vector=np.array([True]*N_HIT+[False]*(Ntot-N_HIT),dtype=bool)
if TARGET_FRACTION>1.0: TARGET_FRACTION=2*BASE_HIT_RATE

# Glide XP Score
setupfile_xp=CTMDDatabase(SETUP_FILE,allow_load=False)
score_data=np.loadtxt(DOCK_SCORES,dtype=str)
xp_scores=score_data[:,1].astype(float)
ligand_names=score_data[:,0].astype(str)
xp_dict=generate_dict(ligand_names,xp_scores)
setupfile_xp.hardcode_scores(xp_dict)

# Boltz2
setupfile_bz=CTMDDatabase(SETUP_FILE,allow_load=False)
score_data=pd.read_csv(BOLTZ2_FILE)

bz_scores=1/(score_data["probability"].to_numpy().astype(float))
ligand_names=score_data["ID"]
bz_dict=generate_dict(ligand_names,bz_scores)
setupfile_bz.hardcode_scores(bz_dict)

# Sampling
sN=dict()
sH=dict()
sN["CTMD"]=[]
sH["CTMD"]=[]
sN["XPScore"]=[]
sH["XPScore"]=[]
sN["Boltz2"]=[]
sH["Boltz2"]=[]

for _ in tqdm.tqdm(range(N_BOOTSTRAP),ncols=100):
    try:
        selected_ligands=np.concatenate((setupfile.sampleByClass(N_HIT,binder=True),setupfile.sampleByClass(binder=False)))
        ligand_scores=np.array([setupfile.getBootstrappedCTMDScore(str(k),commit_frames=COMMIT_FRAMES,cutoff=CUTOFF,timestep=FRAME_TIMESTEP,frame_limit=FRAME_LIMIT)[0] for k in selected_ligands],dtype=float)
        final_hit_vector,final_ctscores,_=sort_from_CTMD(ligand_scores[:,0],hit_vector,ligand_scores[:,1])

        # CTMD
        n_screened,hit_rates=get_enrichment_factor_trend(final_hit_vector,final_ctscores,preordered=True)
        sN["CTMD"].append(n_screened)
        sH["CTMD"].append(hit_rates)

        # XPScore
        ligand_scores=np.array([setupfile_xp.getBootstrappedCTMDScore(str(k))[0] for k in selected_ligands],dtype=float)
        n_xp,xp_hit_rates=get_enrichment_factor_trend(np.array([setupfile_xp.getClass(k) for k in selected_ligands],dtype=bool),-ligand_scores[:,0])
        sN["XPScore"].append(n_xp)
        sH["XPScore"].append(xp_hit_rates)
        
        # Boltz 2
        ligand_scores=np.array([setupfile_bz.getBootstrappedCTMDScore(str(k))[0] for k in selected_ligands],dtype=float)
        n_bz,bz_hit_rates=get_enrichment_factor_trend(np.array([setupfile_bz.getClass(k) for k in selected_ligands],dtype=bool),-ligand_scores[:,0])
        sN["Boltz2"].append(n_bz)
        sH["Boltz2"].append(bz_hit_rates)

    except KeyboardInterrupt: break

sN["CTMD"]=np.array(sN["CTMD"])
sH["CTMD"]=np.array(sH["CTMD"])
sN["XPScore"]=np.array(sN["XPScore"])
sH["XPScore"]=np.array(sH["XPScore"])

# Plot everything
plt.figure(figsize=(8,6))
plt.title(f"Hit rate trend ({N_HIT} binders, {Ntot-N_HIT} non-binders)",fontsize=24)

# CTMD
xdata=np.mean(sN["CTMD"],axis=0)
xshift=(xdata[1]-xdata[0])/8
y_mean=np.mean(sH["CTMD"],axis=0)
y_err=np.std(sH["CTMD"],axis=0)
plt.errorbar(xdata,(y_mean/BASE_HIT_RATE),yerr=(y_err/BASE_HIT_RATE),color="purple",linewidth=3,label="CTMD")

ilocx=np.argmin(np.abs(xdata-TARGET_FRACTION))
print(f"CTMD Enrichment factor at screened fraction {TARGET_FRACTION}:",xdata[ilocx],y_mean[ilocx]/BASE_HIT_RATE,"+/-",y_err[ilocx]/BASE_HIT_RATE,flush=True)

# XPScore
xdata=np.mean(sN["XPScore"],axis=0)
y_mean=np.mean(sH["XPScore"],axis=0)
y_err=np.std(sH["XPScore"],axis=0)
plt.errorbar(xdata+xshift*2,(y_mean/BASE_HIT_RATE),yerr=(y_err/BASE_HIT_RATE),color="blue",linewidth=3,label=DOCK_LABEL)

ilocx=np.argmin(np.abs(xdata-TARGET_FRACTION))
print(f"DOCK Enrichment factor at screened fraction {TARGET_FRACTION}:",xdata[ilocx],y_mean[ilocx]/BASE_HIT_RATE,"+/-",y_err[ilocx]/BASE_HIT_RATE,flush=True)

# Boltz 2
xdata=np.mean(sN["Boltz2"],axis=0)
y_mean=np.mean(sH["Boltz2"],axis=0)
y_err=np.std(sH["Boltz2"],axis=0)
plt.errorbar(xdata+xshift*2,(y_mean/BASE_HIT_RATE),yerr=(y_err/BASE_HIT_RATE),color="green",linewidth=3,label="Boltz2")

ilocx=np.argmin(np.abs(xdata-TARGET_FRACTION))
print(f"Boltz2 Enrichment factor at screened fraction {TARGET_FRACTION}:",xdata[ilocx],y_mean[ilocx]/BASE_HIT_RATE,"+/-",y_err[ilocx]/BASE_HIT_RATE,flush=True)
print()

plt.axhline(y=1.0,linestyle="--",linewidth=3,alpha=0.6,label="Baseline",color='k')
plt.ylabel("Enrichment Factor",fontsize=21)
plt.xlabel("Fraction picked",fontsize=21)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=21)
plt.xlim(TOPK_MIN/Ntot,0.95)
plt.ylim(-0.1,4.5)
plt.show()

# Plot AUROC

# Plot everything
plt.figure(figsize=(8,6))
plt.title(f"Hit rate trend ({N_HIT} binders, {Ntot-N_HIT} non-binders)",fontsize=24)

# CTMD
xdata=np.mean(sN["CTMD"],axis=0)
xshift=(xdata[1]-xdata[0])/8
y_mean=np.mean(sH["CTMD"],axis=0)
y_err=np.std(sH["CTMD"],axis=0)
alphas_CTMD,bedroc_CTMD=getBEDROCProfile(y_mean)
plt.errorbar([0.]+list(xdata),[0.]+list((y_mean*np.arange(1,len(xdata)+1))/N_HIT),yerr=[0.]+list((y_err*np.arange(1,len(xdata)+1))/N_HIT),color="purple",linewidth=3,label="CTMD")

# Glide XP
xdata=np.mean(sN["XPScore"],axis=0)
y_mean=np.mean(sH["XPScore"],axis=0)
y_err=np.std(sH["XPScore"],axis=0)
alphas_xp,bedroc_xp=getBEDROCProfile(y_mean)
plt.errorbar([0.]+list(xdata),[0.]+list((y_mean*np.arange(1,len(xdata)+1))/N_HIT),yerr=[0.]+list((y_err*np.arange(1,len(xdata)+1))/N_HIT),color="blue",linewidth=3,label=DOCK_LABEL)

# Boltz 2
xdata=np.mean(sN["Boltz2"],axis=0)
y_mean=np.mean(sH["Boltz2"],axis=0)
y_err=np.std(sH["Boltz2"],axis=0)
alphas_bz,bedroc_bz=getBEDROCProfile(y_mean)
plt.errorbar([0.]+list(xdata),[0.]+list((y_mean*np.arange(1,len(xdata)+1))/N_HIT),yerr=[0.]+list((y_err*np.arange(1,len(xdata)+1))/N_HIT),color="green",linewidth=3,label="Boltz 2")

plt.plot([0.,1.],[0.,1.],color="grey",linestyle="--",label="Random")
plt.ylabel("Fraction hits found",fontsize=21)
plt.xlabel("Fraction picked",fontsize=21)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=21)
plt.xlim(0.01,0.99)
#plt.ylim(-0.1,3.25)
plt.show()

# Plot BEDROC
plt.figure(figsize=(8,6))
plt.title(f"BEDROC ({N_HIT} binders, {Ntot-N_HIT} non-binders)",fontsize=24)
plt.plot(alphas_CTMD,bedroc_CTMD,color="purple",linewidth=3,label="CTMD")
plt.plot(alphas_xp,bedroc_xp,color="blue",linewidth=3,label=DOCK_LABEL)
plt.plot(alphas_bz,bedroc_bz,color="green",linewidth=3,label="Boltz 2")
plt.legend(fontsize=21)
plt.ylabel("BEDROC Score",fontsize=21)
plt.xlabel("BEDROC $\\alpha$",fontsize=21)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

np.save(BASE_DIR+"/BEDROC_data.npy",(alphas_CTMD,bedroc_CTMD,bedroc_xp,bedroc_bz)) # CTMD, DOCK, Boltz-2
