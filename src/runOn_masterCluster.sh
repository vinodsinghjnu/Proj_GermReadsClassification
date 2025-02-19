#!/bin/bash
#
#SBATCH --ntasks-per-node=1
#SBATCH -J NoJobName
#SBATCH --output=xx-%x.%j.out
#SBATCH --error=xx-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xx@gmail.com
#SBATCH --verbose


## Shift first argument (script path)

#echo "Script path: $script_path"
#echo "Arguments: $@"

#script_path=$1
conda_env=$1
script_path=$2
shift

echo "======"
# Determine script type based on file extension
file_extension="${script_path##*.}"

# Run script based on extension, passing remaining arguments
case "$file_extension" in
    "R")
        source /home/vinodsingh/miniforge3/bin/activate $conda_env
        echo "Cond env1: $CONDA_DEFAULT_ENV"
        which R
        #Rscript --vanilla "$script_path" $@
        Rscript --vanilla $@
        ;;
    "py")
        source /home/vinodsingh/miniforge3/bin/activate $conda_env
        echo "Cond env1: $CONDA_DEFAULT_ENV"
        which python
        #echo  python $@
        python -u $@
        ;;
    "sh")
        source /home/vinodsingh/miniforge3/bin/activate $conda_env
        echo "Cond env1: $CONDA_DEFAULT_ENV"
        echo "COMMAND SUBMITTED"
        echo sh $@
        sh $@
        ;;
    *)
        echo "Unsupported script type: $file_extension"
        exit 1
        ;;
esac


# Print out values of the current jobs SLURM environment variables
env | grep SLURM
# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch


