rule index_picard:
    input:
        ref="data/references/{run}/{sample}_ref.fa",
        checkpoint="data/references/{run}/{sample}_ref.fa.prep",
    output:
        ref_dict="data/references/{run}/{sample}_ref.dict",
    conda:
        "../envs/snp.yaml"
    log:
        "data/references/{run}/logs/{sample}.reference-picard.log",
    resources:
        mem_mb=2000,
#    cache: True
    shell:
        '''
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        echo $PICARD_BIN

        java -Xmx${{mem_gb}}G -jar $PICARD_BIN CreateSequenceDictionary R={input.ref} O={output.ref_dict} &> {log}
        '''

    # run:
    #     if hasattr(resources, "mem_mb"):
    #         mem_gb = resources.mem_mb // 1000
    #     else:
    #         mem_gb = 2

    #     shell(

rule index_fai:
    input:
        ref="data/references/{run}/{sample}_ref.fa",
        checkpoint="data/references/{run}/{sample}_ref.fa.prep",
    output:
        "data/references/{run}/{sample}_ref.fa.fai",
    conda:
        "../envs/snp.yaml"
    shell:
        "samtools faidx {input.ref}"


rule index_bwa:
    input:
        ref="data/references/{run}/{sample}_ref.fa",
        checkpoint="data/references/{run}/{sample}_ref.fa.prep",
    output:
        multiext("data/references/{run}/{sample}_ref.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    conda:
        "../envs/snp.yaml"
    shell:
        "bwa index {input.ref}"


rule index_minimap2:
    input:
        ref="data/references/{run}/{sample}_ref.fa",
        checkpoint="data/references/{run}/{sample}_ref.fa.prep",
    output:
        "data/references/{run}/{sample}_ref.fa.mmi",
    conda:
        "../envs/assembly.yaml"
    shell:
        "minimap2 -d {output} {input.ref}"

