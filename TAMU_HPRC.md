## TAMU-HPRC Job Script Structure
~~~
#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc_maize           # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file
## These are the parameters for TAMU HPRC slurm script.
~~~
### Job submission and monitoring job
* Submit a job: **sbatch** [script_name]
* Cancel/Kill a job: **scancel** [job_id]
* Check status of a single job: **squeue** --job [job_id]
### Searching Software
To search a software in the grace cluster, this [HPRC Available Software webpage](https://hprc.tamu.edu/kb/Software/) needs to be visited and type the software name. The following codes need to be used often to search, load, and unload the software in the HPRC terminal:
* Finding most available software on Grace: **module avail**
* Searching particular software (if I know the software name): **module spider**
* Loading the software: **module load** *software name*
* To know how many softwares or modules are loaded in the current terminal: **module list**
* To remove loaded modules: **module purge**
> It is always recommended to use module purge before using another modules in the same terminal session 
