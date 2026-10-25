#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=gapit_cmlm      # job name
#SBATCH --time=72:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=24           # CPUs (threads) per command
#SBATCH --mem=340G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout_cmlm.%j       # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr_cmlm.%j        # save stderr to file

# load the modules
ml purge
ml WebProxy
ml GCC/13.3.0 OpenMPI/5.0.3 R_tamu/4.4.2 R-bundle-CRAN/2024.11 R-bundle-Bioconductor/3.20-R-4.4.2

# Set the R library path 
# No longer needed since --rtamuenvs flag will take care of this
#export R_LIBS_USER="/scratch/user/farid-bge/R_LIBS_USER/myRlibs"

# run the R script
Rscript --rtamuenvs=GAPIT script_GAPIT_cmlm.R

echo "Run Finish_cmlm"





