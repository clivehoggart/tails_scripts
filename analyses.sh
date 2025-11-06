#BSUB -L /bin/sh
#BSUB -n 3
#BSUB -J CandT[1-4600]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -W 1:00
#BSUB -P acc_psychgen
#BSUB -o CandT.o%J.%I
#BSUB -eo CandT.e%J.%I
#BSUB -M 40000

module load R

Rscript --vanilla /hpc/users/hoggac01/tails_scripts/analyse.R ${LSB_JOBINDEX}
