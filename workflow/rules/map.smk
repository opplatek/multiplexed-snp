rule map_bwamem:
    input:
        r1="data/prep/{run}/{sample}_L001_R1_001_prep.fastq.gz",
        r2="data/prep/{run}/{sample}_L001_R2_001_prep.fastq.gz",
        ref="data/references/{run}/{sample}_ref.fa",
        ref_index=multiext("data/references/{run}/{sample}_ref.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
#        bam=temp(OUTDIR + "{run}/bams/{sample}.bam"),
        bam=OUTDIR + "{run}/bams/{sample}.bam",
        flagstat=OUTDIR + "{run}/qc/mapping_statistics/{sample}-flagstat.txt",
    conda:
        "../envs/snp.yaml"
    threads: 32
    resources:
        tmpdir=lambda wildcards: OUTDIR + f"{wildcards.run}/bams/{wildcards.sample}",
    params:
#        instrument="Illumina"
        extra=r"-R '@RG\tID:{sample}\tLB:{sample}\tSM:{sample}\tPU:{sample}\tPL:Illumina'",
    shell:
        '''
        echo {resources.tmpdir}

        bwa mem -t {threads} -T 20 -v 1 -L 100,100 -M {params.extra} {input.ref} {input.r1} {input.r2} \
        | samtools view -@ {threads} -b - | samtools sort -n -@ {threads} -T {resources.tmpdir} - \
        | samtools fixmate -O bam - - | samtools sort -@ {threads} - > {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        '''
