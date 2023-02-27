# WORK IN PROGRESS!!!
# WORK IN PROGRESS!!!
# WORK IN PROGRESS!!!

# Snakemake workflow: `multiplexed-snp`
A Snakemake workflow for *variant calling for runs with multi-sample/multi-target indexes*

## What does it do
* Calls variants (SNPs and short indels) in samples with multiple genes/genomic regions (for example, amplicons, but also host-parasite samples) multiplexed into a single sequencing index. Carefully aligns reads to the reference, removes overly soft-clipped reads, optionally removes reads with too many mismatches, and calls SNPs/indels. The goal is to minimize false-positive results for the price of missing some true-positive results. 
* The workflow also runs *de novo* assembly and compares the scaffolds with the reference to identify potentially highly mutated samples.
* This workflow is useful for any project where you analyze multiple pooled sequences, which might be somehow similar and you need to be sure there is no *cross-mapping* between samples or references.

### Disclaimer
This is an old(er) pipeline. Some parts are likely outdated and deserve to be updated or completely replaced. However, our collaborators prefer *consistent* results processed the same way as this is a long-term project. Overall, the results we obtained from the workflow seem to be correct.

### Disclaimer 2
This is a highly *experimental* Snakemake workflow. Some parts might be redundant, repetitive, or inefficient. Feel free to make suggestions on how to improve the code or how to make it nicer.

## The workflow
The workflow consists of the following steps:

### Initial QC

### Read preprocessing

### *De novo* assembly

### Mapping

### Mapping postprocessing

### Post-alignment QC

### Variant calling

### Read-based phasing/Haplotype assembly

### Variant summarization

## How to use the workflow
A brief manual on how to use the workflow.

### Environment
#### Primary environment
You have to have [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)/[Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed. Our recommendation is [**Miniconda**](https://docs.conda.io/en/latest/miniconda.html) with [**Mamba**](https://mamba.readthedocs.io/en/latest/installation.html).

With Conda/Miniconda installed, install the primary environment:

```
mamba env create -f environment.yaml -n multiplexed-snp
conda activate multiplexed-snp
```

The absolute minimum you need to run the workflow is Python3 (tested on version 3.7.1) and [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (tested on version >=6.5.2). The rest of the environment in `environment.yaml` are Python packages Snakemake uses to create reports. These are not crucial but help.

#### Rules environments
We recommend to pre-install rules environments before the analysis. This makes it easier to catch installation errors. 

```
threads=6

snakemake --use-conda --conda-frontend "mamba" --conda-create-envs-only -c $threads -p
```

##### Fix `ngsutils`
The latest `ngsutils` release contains an error. Because the tool is not maintained anymore, we have to fix it:

```
# Find where Snakemake/Conda installed ngsutils
filter_fix=$(find .snakemake/conda -name 'filter.py' -print0 | grep -FzZ 'ngsutils/bam/filter.py') # Might give warning "-bash: warning: command substitution: ignored null byte in input"

# Make the fix
echo "Fixing $filter_fix"
sed -i 's/return "maximum mismatch ratio: %s" % self.val/return "maximum mismatch ratio: %s" % self.ratio/' $filter_fix # Fix issue in ngsutils 0.5.9 on line 710 (you can grep -n "return \"maximum mismatch ratio: %s\" % self.val" to be sure); git commit a7f08f5
```

### Inputs

#### Run
The run has to be specified by its **run ID** in `config/runs.csv`. At this moment (02/24/2023), it is not possible to combine multiple runs.

#### Samples
**Raw FASTQ** files must be copied to `data/raw/runID` directory where `runID` is the actual run ID. At this moment (02/24/2023), the workflow can only work with **paired-end** reads.
Samples' names need to be in `sampleID_L001_R1_001.fastq.gz` and `sampleID_L001_R2_001.fastq.gz` format, where `sampleID` is the actual sample ID.

Only samples in `config/samples.csv` are going to be processed. Please edit the config file accordingly.

The workflow has been primarily tested on reads *150 nt* or longer. If you have shorter reads, please decrease the minimum read length (`min_len_fastq`) specified in `config/config.yaml`.

#### References
References must be copied to `data/references/rundID` where `runID` is the same for the samples you want to process. References' names need to be in `sampleID.fa` where `sampleID` is the same sample ID as for the reads.

### Workflow settings
Most of the primary settings can be found in config files defined in `config/config.yaml`. 

#### Mapping filters
The primary filtering settings are:
* Maximum **ratio** of **soft-clipped** relative to the mapped bases (default: 0.1 of soft-clipped bases)
* Maximum **number** of mismatches (default: 999 - no maximum)
* Maximum **ratio** of mismatches relative to the mapped bases (default: 0.1 of mismatched bases)

The settings can be adjusted by editing `config/map_filters.csv`.

#### *De-novo* assembly
The primary settings are:
* Minimum length of the assembled scaffolds (default: 300 nt). You can use `1` to get all lengths.
* Number of the top*N* longest scaffolds to extract (default: 30). You can use `All` to get all scaffolds

The settings can be adjusted by editing `config/assembly.csv`.

### Run
When both references and samples are ready, you can run the workflow:

```
threads=6
export TMPDIR="." # Temp dir Snakemake will use
date=$(date +"%Y%d%d_%H%M%S")

snakemake --use-conda \
    -c $threads --resources assembly_jobs=1 \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.multiplexed-snp.snakemake-stats.txt \
    -p
```

The `--resources assembly_jobs=1` limits the number of concurrently *de novo* assembly jobs to one. *De novo* assembly tends to be RAM-hungry even on small(er) amplicons. You can also add `--set-threads assembly_spades=$threads` to force every *de novo* assembly job to use all the threads, limiting it to a maximum of one job at a time.

You can add `-n` to the main Snakemake command to make a *dry-run* to see what commands would be executed without actually running the workflow.


Optionally, you can create reports and summaries:

```
# Report
snakemake \
    -c $threads \
    --report workflow/report/${date}.multiplexed-snp.snakemake-report.html 

# Summary
snakemake \
    -c $threads \
    --detailed-summary > workflow/report/${date}.multiplexed-snp.snakemake-summary.txt
```

If you need to rerun the whole analysis, you have to remove the outputs first and then rerun Snakemake:

```
threads=6

snakemake --delete-all-output -c $threads
```

You can also run only incomplete jobs in case some of the jobs did not finish:

```
threads=6
date=$(date +"%Y%d%d_%H%M%S")

snakemake --use-conda \
    -c $threads --resources assembly_jobs=1 \
    --configfile config/config.yaml \
    --stats workflow/report/${date}.multiplexed-snp.snakemake-stats.txt \
    -p \
    --rerun-incomplete
 ```

### Output
The default results directory is `results/runID`. 

#### Output directories
The output files are structured:
* **assembly** - **_de novo_ assembly** results
    * **sampleID** - individual sample's assemblies
    * **logs** - log files
    * **scaffolds** - filtered scaffolds (min. length; top*N* longest)
        * **bams** - scaffolds alignment to the reference
* **bams** - **mapping** results **before** any filtering
* **qc** - **QC** results
* **sfPct10mmNo999mmPct10** - this is an example of **variant call** results with a maximum of *10%* soft-clipped bases, a maximum of *999* mismatches (number), and a maximum of *10%* mismatches
    * **bams** - **mapping** results **after** filtering (max. ratio soft-clipping, max. number of mismatches, max. ratio of mismatches)
    * **variants** - variant call results
        * **sampleID** - individual sample's variant calls
        * **extract** - summary/overlap table of variant calls of the different variant callers
        * **gatk** - GATK results
        * **gt** - summary/overlap table of genotyping results of the different variant callers
        * **samtools** - Samtools results
        * **snp** - positions and read counts at variant sites
        * **varscan2** - VarScan2 results
        * **vcf** - VCF file of the variant calls

#### Main output files

#### Supplementary output files