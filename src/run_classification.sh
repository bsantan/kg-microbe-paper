#!/bin/bash
#SBATCH --job-name=classification_analysis
#SBATCH --account=m4689
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256GB
#SBATCH --time=04:00:00
#SBATCH --output=Classification_gold_standard_comparison_%j.out
#SBATCH --error=Classification_gold_standard_comparison_%j.err

# Load required modules
module load python

# Activate virtual environment
source /global/cfs/cdirs/m4689/brook/kg-microbe-paper/src/venv_projects/bin/activate

# Navigate to working directory
cd /global/cfs/cdirs/m4689/kg-microbe-paper/src

# Run the script
echo "Starting Classification_gold_standard_comparison.py at $(date)"
echo "Running on node: $HOSTNAME"
echo "Memory allocated: 256GB"
echo ""

time python Classification_gold_standard_comparison.py

echo ""
echo "Completed at $(date)"
