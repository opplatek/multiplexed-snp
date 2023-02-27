#!/bin/bash
#
# Example of installation and running of Multiplexed SNP Snakemake worfklow
#

### Set default parameters
threads=6

date=$(date +"%Y%d%d_%H%M%S")

export TMPDIR="." # Temp dir Snakemake will use

### Install environments - Run only once!
# source ~/tools/miniconda3/bin/activate # Temporary
mamba env create -f environment.yaml -n multiplexed-snp

conda activate multiplexed-snp
conda config --add channels defaults
conda config --add channels r
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install and fix ngsutils
# TODO: Find a different way how to do this so we can skip ngsutils completely; Will be preinstalled in docker
snakemake --use-conda --conda-frontend "mamba" --conda-create-envs-only -c $threads -p # Preinstall environments including ngsutils env
filter_fix=$(find .snakemake/conda -name 'filter.py' -print0 | grep -FzZ 'ngsutils/bam/filter.py') # Might give warning "-bash: warning: command substitution: ignored null byte in input"
echo "Fixing $filter_fix"
sed -i 's/return "maximum mismatch ratio: %s" % self.val/return "maximum mismatch ratio: %s" % self.ratio/' $filter_fix # Fix issue in ngsutils 0.5.9 on line 710 (you can grep -n "return \"maximum mismatch ratio: %s\" % self.val" to be sure); git commit a7f08f5

### Check the workflow "grammar"
# snakemake --lint # Check correct Snakemake "grammar"

### Run the workflow
## Note: We can "force" assembly (in this case) to prioritize no. of cores for one job over multiple jobs with less cores
# snakemake --use-conda -c $threads --set-threads assembly_spades=$threads -p

# snakemake --delete-all-output -c $threads # Delete all outputs if we want to completely restart the workflow

# Workflow
snakemake --use-conda \
    -c $threads --resources assembly_jobs=1 \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.multiplexed-snp.snakemake-stats.txt \
    -p
# --rerun-incomplete # Rerun incomplete stages

# Html report
snakemake \
    -c $threads \
    --report workflow/report/${date}.multiplexed-snp.snakemake-report.html 

# Text summary
snakemake \
    -c $threads \
    --detailed-summary > workflow/report/${date}.multiplexed-snp.snakemake-summary.txt

