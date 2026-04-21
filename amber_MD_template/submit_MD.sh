#!/bin/bash
#SBATCH --job-name=ET1flC_8xvh
#SBATCH --output=debugMetaD.out
#SBATCH --error=debugMetaD.err
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32000
#SBATCH --account=tiwary-prj-paid


ACCT="tiwary-prj-paid"
export SLURM_ACCOUNT=$ACCT
export SBATCH_ACCOUNT=$ACCT
export SALLOC_ACCOUNT=$ACCT

#source /home/venkata/group_shared/compiled_modules/ambertools24_plmd/sourceme.sh
export AMBERTOOLS24_IMAGE_PATH=/home/venkata/singularity_images/ambertools24.sif
folname=`basename $(pwd)`
echo $folname
singularity exec --nv --bind ../:/scratch $AMBERTOOLS24_IMAGE_PATH /scratch/$folname/amber_MD_template/singularity_passer.sh $folname 0

# Setup reference for PLUMED
cpptraj -i amber_MD_template/extract.in
vmd -dispdev text -e amber_MD_template/weightref.tcl
cp amber_MD_template/plumed_template.dat plumed_bpmetad.dat # Change this line to change how plumed is set up

# Actual BPMetad
singularity exec --nv --bind ../:/scratch $AMBERTOOLS24_IMAGE_PATH /scratch/$folname/amber_MD_template/singularity_passer.sh $folname 1
