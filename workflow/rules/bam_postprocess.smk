rule bam_filter:
    input:
        OUTDIR + "{run}/bams/{sample}.bam",
    output:
        sorted=temp([OUTDIR + "{run}/bams/{sample}.sorted.bam", OUTDIR + "{run}/bams/{sample}.sorted.bam.bai"]),
        filt_sc=temp([OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt-sc.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt-sc.bam.bai"]),
        filt_minlen=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt-minlen.bam"),
        filt_mm=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt-mm.bam"),
        filt=temp([OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt.bam.bai"]),
    conda:
        "../envs/ngsutils_rseqc.yaml"
    params:
        min_len=150,
        max_len=300,
        min_len_map=70,
#        max_softlip_perc=0.1,
#        max_mism_number=999,
#        max_mism_perc=0.1,
        max_softlip_perc=MAP_FILTER_PARAMS["max_softlip_perc"],
        max_mism_number=int(MAP_FILTER_PARAMS["max_mism_number"]),
        max_mism_perc=MAP_FILTER_PARAMS["max_mism_perc"],
    resources:
        mem_mb=2000,
    log:
        OUTDIR + "{run}/qc/" + OUTDIR_SUB + "bam_filter/{sample}.log",
    shell:
        '''
        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`
        echo $mem_gb

        samtools view -h -@ {threads} -F 12 -f 2 -F 2304 -q 10 {input} \
            | grep -v "SA:" \
            | samtools view -@ {threads} -b - \
            | java -Xmx${{mem_gb}}G -jar $PICARD_BIN SortSam SORT_ORDER=coordinate INPUT=/dev/stdin OUTPUT={output.sorted[0]}
        samtools index {output.sorted[0]}

        tmp_file=$(realpath "{input}")

        bamutils removeclipping {output.sorted[0]} - \
            | samtools view -@ {threads} - \
            | grep 'ZC:f:' | awk '{{for (i=1;i<=NF;i++){{if ($i ~/ZC:f:/) {{print $1, $i}}}}}}' | sed 's/ZC:f://' \
            | awk -v MAX_SOFTCLIP="{params.max_softlip_perc}" ' $2 > MAX_SOFTCLIP {{print $1}}' \
            | cut -f1 | awk '!x[$0]++' > ${{tmp_file}}.toremove_highsoftclip.txt

        if [ -s ${{tmp_file}}.toremove_highsoftclip.txt ]; then
            echo "There is `wc -l ${{tmp_file}}.toremove_highsoftclip.txt | cut -d' ' -f1` too much soft-clipped reads with {params.max_softlip_perc} soft-clipping"
            java -jar $PICARD_BIN FilterSamReads I={output.sorted[0]} O={output.filt_sc[0]} READ_LIST_FILE=${{tmp_file}}.toremove_highsoftclip.txt FILTER=excludeReadList
        else
            echo "There are none too much soft-clipped reads with {params.max_softlip_perc} soft-clipping. Continue without filtering."
        fi > {log}
        rm ${{tmp_file}}.toremove_highsoftclip.txt
        samtools index {output.filt_sc[0]}

        reformat.sh in={output.filt_sc[0]} out={output.filt_minlen} minlength={params.min_len} maxlength={params.max_len} overwrite=true
        bamutils filter {output.filt_minlen} {output.filt_mm} -failed ${{tmp_file}}.removed_mismatch.txt -mapped -nosecondary -mismatch {params.max_mism_number} -maximum_mismatch_ratio {params.max_mism_perc}
        rm ${{tmp_file}}.removed_mismatch.txt

        samtools view -h -@ {threads} {output.filt_mm} \
            | perl -slane '$l = 0; $F[5] =~ s/(\d+)[MX=DN]/$l+=$1/eg; print if $l > $MIN_LENGTH_MAPPED or /^@/' -- -MIN_LENGTH_MAPPED={params.min_len_map} \
            | samtools view -@ {threads} -b - > {output.filt[0]}

        samtools index {output.filt[0]}
        '''

rule bam_realign:
    input:
        bam=[OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt.bam.bai"],
        ref="data/references/{run}/{sample}_ref.fa",
        ref_fai="data/references/{run}/{sample}_ref.fa.fai",
        ref_dict="data/references/{run}/{sample}_ref.dict",
    output:
        interval=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.filt.bam.intervals"),
        realign=[OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.indelRealigned.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.indelRealigned.bam.bai"],
    conda:
        "../envs/snp.yaml"
    resources:
        mem_mb=2000,
        tmpdir=lambda wildcards: OUTDIR + f"{wildcards.run}/" + OUTDIR_SUB + f"bams/{wildcards.sample}",
    shell:
        '''
        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        GATK_BIN=`find $CONDA_PREFIX -type f -name 'GenomeAnalysisTK.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T RealignerTargetCreator --num_threads {threads} -R {input.ref} -I {input.bam[0]} -o {output.interval}
        java -Xmx${{mem_gb}}G -Djava.io.tmpdir={resources.tmpdir} -jar $GATK_BIN -T IndelRealigner -R {input.ref} -I {input.bam[0]} --disable_bam_indexing \
            -LOD 0.4 --consensusDeterminationModel USE_SW -targetIntervals {output.interval} -o {output.realign[0]}
        java -Xmx${{mem_gb}}G -jar $PICARD_BIN BuildBamIndex INPUT={output.realign[0]} OUTPUT={output.realign[1]}
        '''

rule bam_duplicates:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.indelRealigned.bam"
    output:
        [OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.duplRem.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.duplRem.bam.bai"],
    conda:
        "../envs/snp.yaml"
    log:
        OUTDIR + "{run}/qc/postprocess/" + OUTDIR_SUB + "{sample}.dupl_stats.log"
    resources:
        mem_mb=2000,
    shell:
        '''
        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        java -Xmx${{mem_gb}}G -jar $PICARD_BIN MarkDuplicates INPUT={input} OUTPUT={output[0]} METRICS_FILE={log} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=true # Removes ALL duplicates, PCR and optical. If you want just to tag them set REMOVE_DUPLICATES=false
        samtools index {output[0]}
        '''
