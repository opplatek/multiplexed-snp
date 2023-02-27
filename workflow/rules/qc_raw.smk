rule fastqc:
    input:
        "data/raw/{run}/{sample}_L001_R{read}_001.fastq.gz",
    output:
        OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001_fastqc.html",
        OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001_fastqc.zip",
    conda:
        "../envs/qc.yaml"
    params:
        odir=OUTDIR + "{run}/qc/fastqc/raw",
    shell:
        "fastqc --quiet -o {params.odir} {input}"


rule fastqc_reads:
    input:
        OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001_fastqc.zip",
    output:
        temp(OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001.total_sequences.txt"),
    params:
        zip=lambda wildcards, input: ".".join(os.path.basename(input[0]).split(".")[:-1]), # remove suffix from the input
        odir=OUTDIR + "{run}/qc/fastqc/raw",
    shell:
        '''
        unzip -q {input} -d {params.odir}
        echo -ne {params.zip}"\\t" > {output}
        grep "Total Sequences" {params.odir}/{params.zip}/fastqc_data.txt \
            | cut -f2 >> {output}
        rm -r {params.odir}/{params.zip}
        '''


rule fastqc_reads_collect:
    input:
        expand(OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001.total_sequences.txt", run=RUNS, sample=SAMPLES, read=READS),
    output:
        OUTDIR + "{run}/qc/fastqc/raw/total_sequences.txt",
    shell:
        "cat {input} > {output}"


rule multiqc:
    input:
        expand(OUTDIR + "{run}/qc/fastqc/raw/{sample}_L001_R{read}_001_fastqc.zip", run=RUNS, sample=SAMPLES, read=READS),
    output:
        OUTDIR + "{run}/qc/fastqc/raw/multiqc_report.html",
    conda:
        "../envs/qc.yaml"
    log:
        OUTDIR + "{run}/qc/fastqc/raw/logs/multiqc.log",
    params:
        odir=OUTDIR + "{run}/qc/fastqc/raw",
    shell:
        "multiqc -f {input} -o {params.odir} -n multiqc_report 2> {log}"

