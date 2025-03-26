#!/bin/bash
#SBATCH --job-name=checkerboard             # Job name
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=10
#SBATCH --partition=public-bigmem,public-cpu,shared-cpu
#SBATCH --time=02:00:00                     # Time limit (hh:mm:ss)
#SBATCH --mem=200GB                         # Memory allocation (in MB)
#SBATCH --output=/home/users/s/sfalcinj/VERTIN/outslurm/%j.out  # Standard output
##SBATCH --mail-type=END,FAIL                # Notifications for job start, end, and failure

# Set a custom MATLAB log directory
export MATLAB_LOG_DIR=/home/users/h/sfalcinj/matlab/java_log

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node(s): $SLURM_NODELIST"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Number of CPUs allocated: $SLURM_CPUS_ON_NODE"

# Load MATLAB module
module load MATLAB/2022a

# Extract the T value from the command line argument
if [ -z "$1" ]; then
    echo "Error: T value not provided. Usage: sbatch run.sh T=<value>"
    exit 1
fi

# Parse T value
T="${1#T=}"
echo "Running with T = $T"

# Set dx_grid based on the value of T
case "$T" in
    1) dx_grid=0.2 ;;
    2) dx_grid=0.2 ;;
    3) dx_grid=0.3 ;;
    4) dx_grid=0.4 ;;
    5) dx_grid=0.5 ;;
    *) 
        echo "Error: Invalid T value. Valid values are 1, 2, 3, 4, 5."
        exit 1
        ;;
esac

# Export T and dx_grid
export T
export dx_grid

# Log the configuration
echo "T = $T"
echo "dx_grid = $dx_grid"

# Run MATLAB script with updated T and dx_grid values
matlab -nodisplay -nosplash -r "run('A_test.m'); run('B_test.m'); exit;"
