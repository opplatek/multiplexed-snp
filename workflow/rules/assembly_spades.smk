import os

def get_mem_mb(wildcards, threads):
    return 6000 * threads

rule assembly_spades:
    input:
        r1="data/prep/{run}/{sample}_L001_R1_001_prep.fastq.gz",
        r2="data/prep/{run}/{sample}_L001_R2_001_prep.fastq.gz",
    output:
        contigs=OUTDIR + "{run}/assembly/{sample}/contigs.fasta",
        scaffolds_full=OUTDIR + "{run}/assembly/{sample}/scaffolds.fasta",
        dir=directory(OUTDIR + "{run}/assembly/{sample}/intermediate_files"),
    conda:
        "../envs/assembly.yaml"
    log:
        OUTDIR + "{run}/assembly/logs/{sample}-spades.log",
    threads: 16
    resources:
#        mem_mb=lambda wildcards, threads: 6000 * threads,
        mem_mb=get_mem_mb,
        assembly_jobs=4,
    script:
        "../wrappers/spades/wrapper.py"

#    wrapper:
#        "file:workflow/wrappers/spades"
    # run:
    #     if hasattr(resources, "mem_mb"):
    #         mem_gb = resources.mem_mb // 1000
    #         memory_requirements = f" --memory {mem_gb}"
    #     else:
    #         memory_requirements = ""

    #     output_dir = os.path.split(output.contigs)[0]

    #     shell("spades.py --careful --threads {threads} {memory_requirements} \
    #     --tmp-dir {output.dir} \
    #     -k 21,33,55,77,99,127 -1 {input.r1} -2 {input.r2} -o {output_dir} &> {log}")


rule assembly_quast:
    input:
        expand(OUTDIR + "{run}/assembly/{sample}/scaffolds.fasta", run=RUNS, sample=SAMPLES)
    output:
        directory(OUTDIR + "{run}/qc/assembly/quast_results"),
    conda:
        "../envs/assembly.yaml"
    shell:
        "quast.py --threads {threads} --output-dir {output} {input}"


rule assembly_filter_link:
    input:
        OUTDIR + "{run}/assembly/{sample}/scaffolds.fasta",
    output:
        OUTDIR + "{run}/assembly/scaffolds/{sample}.scaffolds.minlen{min_len}.fasta",
    conda:
        "../envs/assembly.yaml"
#    params:
#        min_len=MIN_LEN_ASSEMBLY if 'MIN_LEN_ASSEMBLY' in globals() else None,
    shell:
        '''
        if [ {wildcards.min_len} == 1 ]; then
            ln -sf ../{wildcards.sample}/$(basename {input}) {output}
        else
            seqtk seq -L {wildcards.min_len} {input} > {output}
        fi
        '''
    # run:
    #     link = os.path.basename(input[0]) # Doesn't work if we have only one input and call it without the [0] or the input name

    #     if wildcards.min_len == 1:
    #         shell("ln -sf ../{wildcards.sample}/{link} {output}")
    #     else:
    #         shell("seqtk seq -L {wildcards.min_len} {input} > {output}")


rule assembly_get_topn:
    input:
        OUTDIR + "{run}/assembly/scaffolds/{sample}.scaffolds.minlen{min_len}.fasta",
    output:
        OUTDIR + "{run}/assembly/scaffolds/{sample}.scaffolds.minlen{min_len}-top{topn}.fasta",
    run:
        link = os.path.basename(input[0]) # Doesn't work if we have only one input and call it without the [0] or the input name

        if wildcards.topn == "All":
            shell("ln -sf {link} {output}")
        else:
            shell("awk \"/^>/ {{n++}} n>{wildcards.topn} {{exit}} {{print}}\" {input} > {output}")


rule assembly_align:
    input:
        assembly=OUTDIR + "{run}/assembly/scaffolds/{sample}.scaffolds.minlen{min_len}.fasta",
        ref=["data/references/{run}/{sample}_ref.fa", "data/references/{run}/{sample}_ref.fa.mmi"]
    output:
        OUTDIR + "{run}/assembly/scaffolds/bams/{sample}.scaffolds.minlen{min_len}.bam",
        OUTDIR + "{run}/assembly/scaffolds/bams/{sample}.scaffolds.minlen{min_len}.bam.bai"
    conda:
        "../envs/assembly.yaml"
    threads: 8
    shell:
        '''
        minimap2 -t {threads} -ax asm5 {input.ref[0]} {input.assembly} \
            | samtools view -@ {threads} -bh - \
            | samtools sort -@ {threads} - > {output[0]}
        samtools index -@ {threads} {output[0]}
        '''

