# Used Rule inheritance

#rule fastqc_prep:
#    input:
#        "data/prep/{run}/{sample}_L001_R{read}_001_prep.fastq.gz",
#    output:
#        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.html",
#        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip",
#    conda:
#        "../envs/qc.yaml"
#    params:
#        odir=OUTDIR + "{run}/qc/fastqc/prep",
#    shell:
#        "fastqc --quiet -o {params.odir} {input}"


use rule fastqc as fastqc_prep with:
    input:
        "data/prep/{run}/{sample}_L001_R{read}_001_prep.fastq.gz",
    output:
        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.html",
        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip",
    params:
        odir=OUTDIR + "{run}/qc/fastqc/prep",


#rule fastqc_prep_reads:
#    input:
#        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip",
#    output:
#        temp(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001.total_sequences.txt"),
#    params:
#        zip=lambda wildcards, input: ".".join(os.path.basename(input[0]).split(".")[:-1]), # remove suffix from the input
#        odir=OUTDIR + "{run}/qc/fastqc/prep",
#    shell:
#        '''
#        unzip -q {input} -d {params.odir}
#        echo -ne {params.zip}"\\t" > {output}
#        grep "Total Sequences" {params.odir}/{params.zip}/fastqc_data.txt \
#            | cut -f2 >> {output}
#        rm -r {params.odir}/{params.zip}
#        '''


use rule fastqc_reads as fastqc_prep_reads with:
    input:
        OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip",
    output:
        temp(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001.total_sequences.txt"),
    params:
        zip=lambda wildcards, input: ".".join(os.path.basename(input[0]).split(".")[:-1]), # remove suffix from the input
        odir=OUTDIR + "{run}/qc/fastqc/prep",


#rule fastqc_prep_reads_collect:
#    input:
#        expand(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001.total_sequences.txt", run=RUNS, sample=SAMPLES, read=READS),
#    output:
#        OUTDIR + "{run}/qc/fastqc/prep/total_sequences.txt",
#    shell:
#        "cat {input} > {output}"


use rule fastqc_reads_collect as fastqc_prep_reads_collect with:
    input:
        expand(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001.total_sequences.txt", run=RUNS, sample=SAMPLES, read=READS),
    output:
        OUTDIR + "{run}/qc/fastqc/prep/total_sequences.txt",


#rule multiqc_prep:
#    input:
#        expand(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip", run=RUNS, sample=SAMPLES, read=READS),
#    output:
#        OUTDIR + "{run}/qc/fastqc/prep/multiqc_report.html",
#    conda:
#        "../envs/qc.yaml"
#    log:
#        OUTDIR + "{run}/qc/fastqc/prep/logs/multiqc.log"
#    params:
#        odir=OUTDIR + "{run}/qc/fastqc/prep",
#    shell:
#        "multiqc -f {input} -o {params.odir} -n multiqc_report 2> {log}"


use rule multiqc as multiqc_prep with:
    input:
        expand(OUTDIR + "{run}/qc/fastqc/prep/{sample}_L001_R{read}_001_prep_fastqc.zip", run=RUNS, sample=SAMPLES, read=READS),
    output:
        OUTDIR + "{run}/qc/fastqc/prep/multiqc_report.html",
    log:
        OUTDIR + "{run}/qc/fastqc/prep/logs/multiqc.log"
    params:
        odir=OUTDIR + "{run}/qc/fastqc/prep",
