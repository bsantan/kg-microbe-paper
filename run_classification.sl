#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --constraint=cpu
#SBATCH --job-name=kg_classification
#SBATCH --output=logs/classification_%j.out
#SBATCH --error=logs/classification_%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Load any required modules (if needed)
# module load python

# Print job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# Run the script
time uv run python src/Classification_gold_standard_comparison.py

# Print completion
echo "Job completed at: $(date)"
