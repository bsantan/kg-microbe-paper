#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --constraint=cpu
#SBATCH --job-name=kg_classification_opt
#SBATCH --output=logs/classification_%j.out
#SBATCH --error=logs/classification_%j.err

# Create logs directory if it doesn't exist
mkdir -p logs

# Load any required modules (if needed)
# module load python

# Set uv cache to local scratch to avoid NFS locking issues
export UV_CACHE_DIR="${TMPDIR:-/tmp}/uv-cache-$$"
export UV_LINK_MODE=copy
mkdir -p "$UV_CACHE_DIR"

# Print job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"
echo "UV cache directory: $UV_CACHE_DIR"

# Run the optimized script (uses DuckDB queries to avoid loading 30GB into pandas)
time uv run python src/Classification_gold_standard_comparison_optimized.py

# Cleanup local cache
rm -rf "$UV_CACHE_DIR"

# Print completion
echo "Job completed at: $(date)"
