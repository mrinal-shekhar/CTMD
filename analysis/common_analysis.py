#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import pandas as pd
from rdkit.ML.Scoring import Scoring

def get_enrichment_factor_trend(hit_vector,scores,preordered=False):
    assert len(hit_vector)==len(scores), "Scores and hit vector sizes don't match"
    if preordered: order=np.arange(len(hit_vector))
    else: order=np.argsort(-scores) # Largest score first

    hit_vector=np.array(hit_vector)[order]
    scores=np.array(scores)[order]
    n_found=0
    hit_rates=[]
    n_screened=[]
    for i in range(len(hit_vector)):
        n_screened.append(i+1)
        if hit_vector[i]: n_found+=1
        hit_rates.append(n_found/(i+1))
    hit_rates=np.array(hit_rates,dtype=float)
    n_screened=np.array(n_screened,dtype=float)/len(hit_vector)
    return n_screened,hit_rates


def processCTMDCOLVAR(cvfile,cutoff=0.35,commit_frames=2,frame_limit=-1):
    #print("Using frame limit:",frame_limit)
    cvdata=np.loadtxt(cvfile,usecols=(0,1,-1))[:frame_limit] # Time, CV, c(t)
    if np.prod(cvdata.shape)<3: return cvdata, -1, np.nan
    if len(cvdata)<commit_frames:
        try: return cvdata, -1, cvdata[-1,-1]
        except: return cvdata, -1, np.nan

    detach=(cvdata[:,1]>cutoff)
    convfilt=np.ones(commit_frames,dtype=float)
    detach=np.convolve(detach,convfilt,mode="valid")
    retidx=np.argmax(detach)
    detqual=(detach[retidx]/commit_frames)
    if detqual<0.99: # Trajectory went to end
        retidx=len(cvdata)-1
    return cvdata,retidx,cvdata[retidx,-1]
def processCTMDFolder(curfol,n_replica=3,cutoff=0.3,commit_frames=200,timestep=0.2,bootstrap_repeats=250,**kwargs): # Timestep is in ps; Take minimum of n_replica repeats. Rest are used as bootstrap basis
    all_colvars=sorted([f for f in os.listdir(curfol) if f.startswith("COLVAR_")])
    if len(all_colvars)<n_replica: raise ValueError("Not enough replicas in folder "+str(curfol))

    cts=[]
    exit_times=[]
    for curcolvar in all_colvars:
        cvdata,exitframe,ct=processCTMDCOLVAR(curfol+"/"+curcolvar,commit_frames=commit_frames,cutoff=cutoff,**kwargs) # 0.2ps/frame = 50ps
        cts.append(ct)
        exit_times.append(timestep*exitframe)

    ctscores=[]
    rtscores=[]
    for _ in range(bootstrap_repeats):
        selcts=np.random.choice(cts,n_replica,replace=False)
        ctscores.append(np.min(selcts))
        selrts=np.random.choice(exit_times,n_replica,replace=False)
        rtscores.append(np.min(selrts))
    ctscore=np.mean(ctscores)
    ctscore_error=np.std(ctscores)
    rtscore=round(np.mean(rtscores),3)
    rtscore_error=np.std(rtscores)
    return ctscore,rtscore,ctscore_error,rtscore_error

def sort_from_CTMD(ctscores,hit_vector,ligand_rtscores=None):
    ligand_ctscores=ctscores
    temp_ligand_order=np.argsort(ligand_ctscores)[::-1] # Higher scores are better for CTMD
    final_ligand_order=[]
    # Resolve ties of <=1kT (2.5 kJ/mol in c(t)) using residence time
    last_ct=10000.0
    last_rt=1e20
    secondary=[]
    nanfails=0
    if ligand_rtscores is None: ligand_rtscores=np.random.normal(0,1,ligand_ctscores.shape)

    for nli in range(len(temp_ligand_order)):
        myct=ligand_ctscores[temp_ligand_order[nli]]
        myrt=ligand_rtscores[temp_ligand_order[nli]]
        if np.isnan(myct):
            nanfails+=1
            myct=0
            myrt=0
        if (last_ct-myct)>2.5: # 1kT = 2.5 kJ/mol
            if len(secondary):
                nord=[ligand_rtscores[i] for i in secondary]
                nord=np.argsort(nord)[::-1]
                for n in nord: final_ligand_order.append(secondary[n])
                secondary=[]
                #print("Secondaries resolved")
            final_ligand_order.append(temp_ligand_order[nli])
            last_ct=myct
            last_rt=myrt
        else:
            if not len(secondary):
                secondary.append(final_ligand_order[-1])
                final_ligand_order=final_ligand_order[:-1]
            secondary.append(temp_ligand_order[nli])
            last_ct=min(myct,last_ct)
            last_rt=min(myrt,last_rt)
    if len(secondary):
        nord=[ligand_rtscores[i] for i in secondary]
        nord=np.argsort(nord)[::-1]
        for n in nord: final_ligand_order.append(secondary[n])

    final_ctscores=[]
    final_hit_vector=[]
    for o in final_ligand_order:
        #print(o,all_folders[o],ligand_ctscores[o],hit_vector[o])
        final_ctscores.append(ligand_ctscores[o])
        final_hit_vector.append(hit_vector[o])
    final_hit_vector=np.array(final_hit_vector,dtype=bool)
    final_ctscores=np.array(final_ctscores,dtype=float)
    return final_hit_vector,final_ctscores,nanfails

class CTMDDatabase:
    def __init__(self,hit_data_file="hit_data_mol2.spc",n_replica=3,allow_load=True):
        self.Nrep=n_replica
        self.data_file=hit_data_file
        data_text=np.loadtxt(self.data_file,dtype=str)
        self.Qnets=data_text[:,2].astype(int)
        self.binder_flag=np.array([(str(v)=="True") or (str(v)!="False" and int(str(v))) for v in data_text[:,1]],dtype=bool) #data_text[:,1].astype(bool)
        self.names=[v.replace(".mol2","") for v in data_text[:,0]]
        self.N=len(self.names)
        self.CTMDdir=os.path.dirname(hit_data_file)
        self.loaded_data=dict()
        self.loadable=allow_load
        print("Detected CTMD directory:",self.CTMDdir,flush=True)

    def getCTMDScore(self,n,return_raw=False,commit_frames=500,cutoff=0.6,timestep=0.2,**kwargs):
        if type(n)==int: n=self.names[n]
        elif type(n)==str: pass
        else: raise ValueError("'n' must be of type int (index) or str (name). Received "+str(type(n)))

        try: i=self.names.index(n)
        except: raise ValueError(f"Name {n} not found in list of ligands loaded!")

        curfol=self.CTMDdir+"/"+n
        if n not in self.loaded_data:
            if not self.loadable: raise ValueError("Data not present for key "+str(n)+" but loading from COLVARs is disabled")
            all_colvars=sorted([f for f in os.listdir(curfol) if f.startswith("COLVAR_")])
            if len(all_colvars)<self.Nrep: raise ValueError("Not enough replicas in folder "+str(curfol))

            cts=[]
            exit_times=[]
            for curcolvar in all_colvars:
                cvdata,exitframe,ct=processCTMDCOLVAR(curfol+"/"+curcolvar,commit_frames=commit_frames,cutoff=cutoff,**kwargs) # 0.2ps/frame = 50ps
                cts.append(ct)
                exit_times.append(timestep*exitframe)
            self.loaded_data[n]=(cts,exit_times)

        cts,exit_times=self.loaded_data[n]
        if return_raw: return cts,exit_times

        selcts=np.random.choice(cts,self.Nrep,replace=False)
        selrts=np.random.choice(exit_times,self.Nrep,replace=False)
        return np.min(selcts),np.min(selrts) # One sample
    def getBootstrappedCTMDScore(self,n,repeats=250,**kwargs):
        stats_raw=np.array([self.getCTMDScore(n,**kwargs) for _ in range(repeats)],dtype=float)
        return np.mean(stats_raw,axis=0),np.std(stats_raw,axis=0)

    def hardcode_scores(self,score_dict,time_dict=None):
        for k in score_dict:
            if k not in self.names: continue
            ct=[score_dict[k]]*self.Nrep
            if time_dict is None: rt=[np.nan]*self.Nrep
            else: rt=[time_dict[k]]*self.Nrep

            self.loaded_data[k]=(ct,rt)

    def sampleByClass(self,n=-1,binder=True,with_replacement=False):
        opts=np.where(self.binder_flag==binder)[0]
        if n<0: idxs=opts
        else: idxs=np.random.choice(opts,n,replace=with_replacement)
        return np.array([self.names[i] for i in idxs],dtype=str)

    def getClass(self,n):
        if type(n)==int: n=self.names[n]
        else: n=str(n)
        return self.binder_flag[self.names.index(n)]
    def forceClass(self,n,cls):
        if type(n)==int: n=self.names[n]
        else: n=str(n)
        self.binder_flag[self.names.index(n)]=bool(cls)
        return bool(cls)
    def eraseLigand(self,n):
        if type(n)==int: n=self.names[n]
        else: n=str(n)
        idx=self.names.index(n)
        self.binder_flag=np.array(list(self.binder_flag[:idx])+list(self.binder_flag[idx+1:]),dtype=bool)
        self.Qnets=np.array(list(self.Qnets[:idx])+list(self.Qnets[idx+1:]),dtype=int)
        self.N-=1
        self.names=self.names[:idx]+self.names[idx+1:]
        if n in self.loaded_data: loaded_data.pop(n)
        return

def getBEDROCProfile(ef_profile):
    ef_profile=np.round(ef_profile*np.arange(1,len(ef_profile)+1))
    hit_no,hit_idx=np.unique(ef_profile,return_index=True)
    if hit_no[0]<0.5: # Initial ligand is non-binder
        hit_no=hit_no[1:]
        hit_idx=hit_idx[1:]

    bitmatch=np.zeros(len(ef_profile),dtype=bool); bitmatch[hit_idx]=True
    v1=np.arange(len(bitmatch))+1
    bedroc_in=np.stack((v1,bitmatch),dtype=float,axis=1)
    bedroc_in=bedroc_in[:500]
    print("BEDROC Scores from",len(bedroc_in),"scores")

    alphas=np.arange(0.1,50,0.1)
    bedroc_scores=[]
    for al in alphas:
        bedroc_scores.append(Scoring.CalcBEDROC(bedroc_in,1,al))
    bedroc_scores=np.array(bedroc_scores)
    print("BEDROC 20:",round(Scoring.CalcBEDROC(bedroc_in,1,20),3))
    return alphas,bedroc_scores
