#!/usr/bin/env python
# coding: utf-8
# author: venkataadury @ github
# From python_libraries (https://github.com/venkataadury/python_libraries/blob/main/python/General.py)
# In[1]:


import numpy as np
from rdkit import Chem


# In[2]:


class SequentialSMILESLoader:
    def __init__(self,smifilename,skip_failures=True,attach_names=False,count_max=True,default_batch=1,max_mols=-1):
        self.filename=smifilename
        if type(smifilename)==str:
            self.file=open(smifilename,"r")
        else:
            self.file=None
        self.counter=0
        self.autoskip=skip_failures
        self.name=attach_names
        self.ended=False
        if count_max:
            if self.file:
                self.linecount=len(self.file.readlines())
                self.file=open(smifilename,"r")
            else:
                self.linecount=len(smifilename) #Assuming input is a list of molecules
        else: self.linecount= None
        self.default_batch=default_batch
        self.mollim=max_mols
        self.skipcount=0
    
    def getFilename(self): return (self.filename if type(self.filename)==str else None)
    def restartSequence(self,keep_location=False):
        ret=SequentialSMILESLoader(self.filename,self.autoskip,self.name,self.linecount is None,self.default_batch)
        if keep_location: ret.skipNext(self.counter)
        return ret
    
    def skipNext(self,num=0):
        K=0
        if self.file is not None:
            for ln in self.file:
                self.counter+=1
                K+=1
                self.skipcount+=1
                if K>=num: break
            else: self.ended=True
        else:
            self.counter+=num
            self.skipcount+=num
            if self.counter>=len(self.filename): self.ended=True
            
    def getNext(self,num=0,as_smiles=False):
        if num==0: num=self.default_batch
        ret=[]
        if self.name:
            names=[]
        K=self.counter
        skip=0
        empty=0
        
        if self.ended: return None
        
        if self.file is None:
            retn=[]
            retm=[]
            while len(retn)<num or num<0:
                self.counter+=1
                if self.counter>=len(self.filename) or (self.mollim>0 and self.counter>=self.mollim):
                    self.ended=True
                    break
                curmol=self.filename[self.counter]
                if type(curmol)==str and (not as_smiles):
                    try:
                        curmol=Chem.MolFromSmiles(curmol)
                    except:
                        continue
                if self.name: retm.append(self.name[self.counter])
                retn.append(curmol)
            if self.name: return retn,retm
            else: return retn
            
        for ln in self.file:
            ln=ln.strip().split()
            if not len(ln):
                self.counter+=1
                empty+=1
                continue
            
            mol=ln[0]
            if self.name: myname=ln[1].strip()
            else: myname = None
            if not as_smiles:
                mol=Chem.MolFromSmiles(mol)
                if self.autoskip and (mol is None):
                    self.counter+=1
                    skip+=1
                    continue
            ret.append(mol)
            if self.name: names.append(myname)
            
            self.counter+=1
            if (self.counter-K-skip-empty==num) or (self.mollim>0 and self.counter-self.skipcount>=self.mollim): break
        else: self.ended=True
        
        if self.name: return ret,names
        else: return ret
    
    def drain(self,as_smiles=False): return self.getNext(-1,as_smiles)
    def __del__(self):
        if self.file is not None:
            self.file.close()
    
    def __help__(self): print("This is help. The sequential smiles loader can (ideally) work with both SMILES files and lists of molecules/smiles strings")


# In[3]:


class RDKitMoleculeWrapper:
    def __init__(self,mol):
        self.mol = mol
        if type(self.mol)!=Chem.rdchem.Mol: self.mol=Chem.MolFromSmiles(self.mol)
    
    def __call__(self): return self.mol
    def get(self): return self.mol
    
    def getSize(self,includeHs=False):
        if includeHs: tempmol=Chem.AddHs(self.mol)
        else: tempmol=self.mol
        return (tempmol.GetNumHeavyAtoms() if (not includeHs) else len(tempmol.GetAtoms()))


# In[4]:


import copy
# Get a similarity score distribution
def getSimilarityDistribution(smilesloader,goal=max,post=None,default=0,fingerprint=Chem.RDKFingerprint,similarity=Chem.DataStructs.TanimotoSimilarity,precompute_fps=False,batch_mols=0):
    num_mol=smilesloader.linecount
    keep_count=False
    if num_mol is None:
        num_mol=0
        keep_count=True
    if type(goal)==str:
        if goal=="max":
            goal=max
            post=None #Post processing
        if goal=="mean":
            goal=lambda x,y:x+y
            post=lambda z: z/num_mol
    
    
    
    if precompute_fps:
        mols=smilesloader.drain(as_smiles=False)
        if type(mols)==tuple and type(mols[0])==list: mols=mols[0]
        mols=[fingerprint(mol) for mol in mols]
        print("Fingerprints computed for",len(mols),"molecules")
    
        ret=[default for _ in mols]
    
        for i in range(len(ret)):
            mol1=mols[i]
            if (i+1)%2500==0: print(i+1,"molecules completed of",len(ret))
            #if not precompute_fps: mol1=fingerprint(mol1)
            for j,mol2 in enumerate(mols):
                if j==i : continue
                #if not precompute_fps: mol2=fingerprint(mol2)
                simsc=similarity(mol1,mol2)
                ret[i]=goal(ret[i],simsc)
    else:
        ret=[]
        mainloader=smilesloader.restartSequence(keep_location=True)
        thismols=mainloader.getNext(batch_mols,as_smiles=False)
        while thismols is not None:
            if type(thismols)==tuple and type(thismols[0])==list: thismols=thismols[0]
            thisfps=[fingerprint(mol) for mol in thismols]
            #Do stuff here
            nxt=[default]*len(thismols)
            shf=len(ret)
            print(shf,"molecules processed")
            startloader=smilesloader.restartSequence(keep_location=True)
            allmols=startloader.getNext(batch_mols,as_smiles=False)
            blockind=0
            while allmols is not None:
                if type(allmols)==tuple and type(allmols[0])==list: allmols=allmols[0]
                allfps=[fingerprint(mol) for mol in allmols]
                for i,tmol in enumerate(thisfps):
                    for j,cfp in enumerate(allfps):
                        if j+blockind==i+shf: continue
                        nxt[i]=goal(nxt[i],similarity(tmol,cfp))
                blockind+=len(allmols)
                allmols=startloader.getNext(batch_mols,as_smiles=False)
            ret+=nxt
            thismols=mainloader.getNext(batch_mols,as_smiles=False)
    
    if post is not None: ret=[post(v) for v in ret]
    return np.array(ret)

def getSimilarityMatrix(smilesloader,default=np.nan,fingerprint=Chem.RDKFingerprint,similarity=Chem.DataStructs.TanimotoSimilarity,precompute_fps=False,batch_mols=0,num_mols=None,force_diagonal=None):
    fullsmiload=copy.deepcopy(smilesloader)
    if num_mols is None:
        if fullsmiload.linecount is None: raise ValueError("Can't have no line counting in the SMILES loader for full matrix. If you want to use this approach, pass no. of molecules as num_mols=... togetSimilarityMatrix")
        else: num_mols=int(smilesloader.linecount)
    if precompute_fps:
        mols=fullsmiload.drain(as_smiles=False) if not fullsmiload.name else fullsmiload.drain(as_smiles=False)[0]
        fps=[Chem.RDKFingerprint(m) for m in mols]
        del mols
        num_mols=len(fps)
        simmat=np.ones((num_mols,num_mols),dtype=float)*default
        for fpi,fp in enumerate(fps):
            simvec=Chem.DataStructs.BulkTanimotoSimilarity(fp,fps)
            assert (simvec[fpi]>=0.999), "Internal error occurred in getSimilarityMatrix. Self-fingerprint of ligand "+str(fpi)+" is not 1"
            if force_diagonal is not None: simvec[fpi]=force_diagonal
            simmat[fpi]=simvec
    else:
        simmat=np.ones((num_mols,num_mols),dtype=float)*default
        print("WARN: Full similarity matrix without precompute_fps will be slow")
        mainloader=smilesloader.restartSequence(keep_location=True)
        
        thismols=mainloader.getNext(batch_mols,as_smiles=False)
        mainindex=0
        while thismols is not None:
            #print(mainindex,"of",num_mols)
            mainfps=[fingerprint(m) for m in thismols]
            subloader=smilesloader.restartSequence(keep_location=True)
            submols=subloader.getNext(batch_mols,as_smiles=False)
            index=0
            while submols is not None:
                subfps=[fingerprint(s) for s in submols]
                for fpi,fp in enumerate(mainfps):
                    simvec=Chem.DataStructs.BulkTanimotoSimilarity(fp,subfps)
                    simmat[mainindex+fpi,index:index+len(subfps)]=simvec
                index+=len(submols)
                submols=subloader.getNext(batch_mols,as_smiles=False)
            mainindex+=len(mainfps)
            thismols=mainloader.getNext(batch_mols,as_smiles=False)
    return simmat

def largest_dissimilar_subset(simmat_,cutoff=0.7,tries=25000,skip=None,silent=False,print_freq=250):
    simmat=simmat_+np.eye(len(simmat_))
    final_sel=[]
    for t in range(tries):
        sel=[]
        allowed=np.ones(len(simmat)).astype(bool)
        if skip is not None: allowed[skip]=False
        idx_list=np.arange(len(simmat))
        while np.any(allowed):
            ch=np.random.choice(idx_list[allowed])
            sel.append(ch)
            #close_idx=np.where()[0]
            allowed[simmat[ch]>=cutoff]=False
        if not silent and (t%print_freq==0): print("Try",t,"picked",len(sel),"ligands")
        if len(sel)>len(final_sel):
            final_sel=sel
    return final_sel


# In[18]:


class Mol2FileLoaded:
    def __init__(self,file,include_comments=True,comment_char='#'):
        if type(file)==str:
            if "\n" in file: file=file.split("\n")
            else: file=list(open(file,"r").readlines())
        elif type(file)==list or type(file)==tuple: pass
        else:
            try: file=list(file.readlines())
            except: print("Assuming input is iterable")

        self.sections=dict()
        self.header_comments=[]
        self.comment_char=comment_char
        self.record_comments=include_comments # All comments move to header
        
        self.load_lines(file)

    def load_lines(self,file):
        sec=None
        for l in file:
            l=l.strip()
            if len(l)<3: continue
            if l[0]==self.comment_char:
                if self.record_comments: self.header_comments.append(l)
                continue
            
            if l.startswith("@<TRIPOS>"):
                l=l.replace("@<TRIPOS>","").strip()
                sec=l
                self.sections[sec]=[]
                continue
            if sec is None: raise ValueError("Non-comment lines before and @<TRIPOS> entries!")
            self.sections[sec].append(l)

    def reconstructMol2Block(self):
        ret=""
        for l in self.header_comments: ret+=l+"\n"
        for s in self.sections:
            ret+="@<TRIPOS>"+s+"\n"
            for l in self.sections[s]: ret+=l+"\n"
        return ret

    def computeNetCharge(self):
        q=0
        for l in self.sections["ATOM"]:
            l=l.strip().split(" ")[-1]
            q+=float(l)
        return q

    def extractName(self): return self.sections["MOLECULE"][0].strip()


# In[51]:


class SequentialMol2Loader:
    def __init__(self,filename,attach_names=False,default_batch=1,max_mols=-1,comment_char='#'):
        self.filename=filename
        if type(filename)==str: self.file=open(filename,"r")
        else:
            raise ValueError("Input to SequentialMol2Loader should be a mol2 filename")
        self.counter=0
        self.name=attach_names
        self.ended=False
        self.default_batch=default_batch
        self.mollim=max_mols
        self.skipcount=0
        self.comment_char=comment_char
        self.new_header=self.skimTop()

    def getFilename(self): return self.filename
    def restartSequence(self,keep_location=False):
        ret=SequentialMol2Loader(self.filename,self.name,self.default_batch)
        if keep_location: ret.skipNext(self.counter)
        return ret

    def skimTop(self):
        ret=[]
        for l in self.file:
            l=l.strip()
            if len(l)<3: continue
            if l[0]==self.comment_char:
                ret.append(l)
                continue
            break
        return ret
    def readNextBlock(self):
        self.header=self.new_header
        self.new_header=[]
        lines=[]
        for l in self.file:
            l=l.strip()
            if len(l)<3: continue
            if l[0]==self.comment_char:
                if len(lines): self.new_header.append(l)
                else: self.header.append(l)
                continue
            if "@<TRIPOS>MOLECULE" in l:
                break
            else: lines.append(l)
        else: self.ended=True
        final_lines=self.header+["@<TRIPOS>MOLECULE"]+lines
        self.counter+=1
        return final_lines

    def skipNext(self,num):
        for _ in range(num): _ = self.readNextBlock()
    def getNext(self,num=0,as_text=False):
        if num==0: num=self.default_batch
        ret=[]
        if self.name: names=[]
        if self.ended: return None
        
        for i in range(num):
            m2b=self.readNextBlock()
            if not as_text: m2b=Mol2FileLoaded(m2b,comment_char=self.comment_char)
            if self.name:
                if as_text: molname=m2b[m2b.find("@<TRIPOS>MOLECULE")+1].strip()
                else: molname=m2b.extractName()
            ret.append(m2b)
            if self.name: names.append(molname)
            if self.mollim>0 and self.counter>=mollim: self.ended=True
        if self.name: return ret,names
        else: return ret


# In[9]:


# Usage for getSimilarityDistribution()
'''
import time

starttime=time.time()
file="/home/venkata/CADD/denovo_benchmarking/analysis/test.smi"
sload=SequentialSMILESLoader(file,default_batch=256)
sload.skipNext(1024) #Skipping Part of the Iterator is Allowed
sdist=getSimilarityDistribution(sload,goal="max",precompute_fps=True) #You can use precompute_fps=False as well to save RAM
endtime=time.time()

import matplotlib.pyplot as plt
plt.hist(sdist)
print(endtime-starttime,len(sdist))
'''
pass


# In[7]:


print("Loaded Library 'General (RDKit)'")

