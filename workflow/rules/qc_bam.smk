rule samtools_flagstat:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.{bamtype}.bam",
    output:
        OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.{bamtype}-flagstat.txt",
    conda:
        "../envs/qc.yaml"
    shell:
        "samtools flagstat {input} > {output}"


rule rseqc:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.{bamtype}.bam",
    output:
        multiext(OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.{bamtype}", ".qual.boxplot.pdf", ".DupRate_plot.pdf", ".GC_plot.pdf", ".NVC_plot.pdf", ".stats"),
    conda:
        "../envs/ngsutils_rseqc.yaml"
    params:
        prefix=OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.{bamtype}",
    shell:
        '''
        read_quality.py -i {input} -o {params.prefix}
        read_duplication.py -i {input} -o {params.prefix}
        read_GC.py -i {input} -o {params.prefix}
        read_NVC.py -i {input} -o {params.prefix}
        bam_stat.py -q 30 -i {input} &> {params.prefix}.stats
        '''


rule qualimap:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.{bamtype}.bam",
    output:
        pdf=OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.{bamtype}-qualimap.pdf",
        genres=OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.{bamtype}-qualimap.genome_results.txt",
    conda:
        "../envs/qc.yaml"
    threads: 12
    shell:
        '''
        unset DISPLAY

        qualimap bamqc -bam {input} -nt {threads} -c -outdir $(dirname {output.pdf})/{wildcards.sample} -outfile $(basename {output.pdf})
        mv $(dirname {output.pdf})/{wildcards.sample}/$(basename {output.pdf}) {output.pdf}
        mv $(dirname {output.pdf})/{wildcards.sample}/genome_results.txt {output.genres}
        '''


rule mappings_get:
    input:
        flagstat_pre=expand(OUTDIR + "{run}/qc/mapping_statistics/{sample}-flagstat.txt", run=RUNS, sample=SAMPLES),
        flagstat_post=expand(OUTDIR + "{run}/qc/" + OUTDIR_SUB + "postprocess/{sample}.duplRem-flagstat.txt", run=RUNS, sample=SAMPLES),
    output:
        map_pre=OUTDIR + "{run}/qc/mapping_statistics/preFilter.flagstat.txt",
        map_post=OUTDIR + "{run}/qc/" + OUTDIR_SUB + "mapping_statistics/postFilter.flagstat.txt",
    shell:
        """
        for flagstat in {input.flagstat_pre}; do
            map_line=`grep \"+ 0 mapped\" $flagstat`
            sample=`echo $(basename $flagstat) | sed 's/-flagstat.txt//'`
            echo -e "$sample\\t$sample\\t$map_line"
        done > {output.map_pre}

        for flagstat in {input.flagstat_post}; do
            map_line=`grep \"+ 0 mapped\" $flagstat`
            sample=`echo $(basename $flagstat) | sed 's/-flagstat.txt//'`
            echo -e "$sample\\t$sample\\t$map_line"
        done > {output.map_post}
        """


rule mappings_plot:
    input:
        OUTDIR + "{run}/qc/" + OUTDIR_SUB + "mapping_statistics/postFilter.flagstat.txt",
    output:
        OUTDIR + "{run}/qc/" + OUTDIR_SUB + "mapping_statistics/postFilter.flagstat.pdf",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/mapping_stats.R"
