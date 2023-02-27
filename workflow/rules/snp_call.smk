rule bam_recalibr:
    input:
        bam=OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.duplRem.bam",
        ref="data/references/{run}/{sample}_ref.fa",
    output:
        bam=[OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.recal.bam", OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.recal.bam.bai"],
        vcf_raw=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.gatk.raw.vcf"),
        vcf_raw_idx=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.gatk.raw.vcf.idx"),
        recal_table=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.table"),
        recal_table_post=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal-post.table"),
        recal_pdf=OUTDIR + "{run}/qc/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.pdf",
    conda:
        "../envs/snp.yaml"
    params:
        call_conf=30.0,
        emit_conf=10.0
    resources:
        mem_mb=2000,
    shell:
        '''
        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        GATK_BIN=`find $CONDA_PREFIX -type f -name 'GenomeAnalysisTK.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        java -Xmx${{mem_gb}}G -jar $PICARD_BIN BuildBamIndex INPUT={input.bam} OUTPUT={input.bam}.bai
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T HaplotypeCaller -I {input.bam} -R {input.ref} -o {output.vcf_raw} -stand_call_conf {params.call_conf} -stand_emit_conf {params.emit_conf}
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {output.vcf_raw} -o {output.recal_table}
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {output.vcf_raw} -BQSR {output.recal_table} -o {output.recal_table_post}
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T AnalyzeCovariates -R {input.ref}  -before {output.recal_table} -after {output.recal_table_post} -plots {output.recal_pdf}
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T PrintReads -I {input.bam} -R {input.ref} -BQSR {output.recal_table} --disable_bam_indexing -o {output.bam[0]}
        java -Xmx${{mem_gb}}G -jar $PICARD_BIN BuildBamIndex INPUT={output.bam[0]} OUTPUT={output.bam[1]}
        '''


rule snp_call_gatk:
    input:
        bam=OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.recal.bam",
        ref="data/references/{run}/{sample}_ref.fa",
    output:
        vcf_raw=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.raw.vcf"),
        vcf_raw_idx=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.raw.vcf.idx"),
        vcf_filt=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.vcf",
        vcf_filt_idx=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.vcf.idx"),
        vcf_filt_phased=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.phased.vcf",
        vcf_filt_phased_idx=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.phased.vcf.idx"),
    conda:
        "../envs/snp.yaml"
    params:
        call_conf=30.0,
        emit_conf=10.0,
        phase_thresh=20.0,
    resources:
        mem_mb=2000,
    shell:
        '''
        PICARD_BIN=`find $CONDA_PREFIX -type f -name 'picard.jar'`
        GATK_BIN=`find $CONDA_PREFIX -type f -name 'GenomeAnalysisTK.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T HaplotypeCaller -I {input.bam} -R {input.ref} -o {output.vcf_raw} \
            -stand_call_conf {params.call_conf} -stand_emit_conf {params.emit_conf}
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T VariantFiltration --variant:VCF {output.vcf_raw} -R {input.ref} -o {output.vcf_filt} \
            --filterExpression "MQ >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" \
            --filterExpression "DP < 5" --filterName "LowCoverage" \
            --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" \
            --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" \
            --filterExpression "QD < 1.5" --filterName "LowQD" --filterExpression "SB > -10.0" --filterName "StrandBias" \
            --filterExpression "FS > 60.0" --filterName "StrandBiasFS"
        java -Xmx${{mem_gb}}G -jar $GATK_BIN -T ReadBackedPhasing -I {input.bam} --variant:VCF {output.vcf_filt} -R {input.ref} -o {output.vcf_filt_phased} \
            --phaseQualityThresh {params.phase_thresh}
        '''


rule snp_call_samtools:
    input:
        bam=OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.recal.bam",
        ref="data/references/{run}/{sample}_ref.fa",
    output:
        pileup=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/samtools/{sample}.recal.mpileup"),
        bcf=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/samtools/{sample}.recal.bcf"),
        vcf_filt=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/samtools/{sample}.recal.samtools.filt.vcf",
    conda:
        "../envs/samtools-0.1.19.yaml"
    params:
        min_mapq=1,
        indel_min_read=3,
        mapq_adj=50,
        max_depth=100000,
    shell:
        '''
        samtools view -b -u -q {params.min_mapq} {input.bam} \
            | samtools mpileup -A -u -B -m {params.indel_min_read} -C {params.mapq_adj} -p -d {params.max_depth} -f {input.ref} - \
            > {output.pileup}
        bcftools view -bvcg {output.pileup} > {output.bcf}
        bcftools view {output.bcf} \
            | vcfutils.pl varFilter -D 10000000 \
            > {output.vcf_filt}
        '''


rule snp_call_varscan2:
    input:
        bam=OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.recal.bam",
        ref="data/references/{run}/{sample}_ref.fa",
    output:
        pileup=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.mpileup2"),
        vcf_snp=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.snp.vcf",
        vcf_indel=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.indel.vcf",
    conda:
        "../envs/samtools-0.1.19.yaml"
    params:
        min_mapq=1,
        indel_min_read=3,
        mapq_adj=50,
        max_depth=100000,
        min_freq=0.05,
        pval=0.05,
        strand=1,
    resources:
        mem_mb=2000,
    shell:
        '''
        VARSCAN2_BIN=`find $CONDA_PREFIX -type f -name 'VarScan.jar'`
        mem_gb=`echo "{resources.mem_mb} / 1000" | bc`

        samtools view -b -u -q {params.min_mapq} {input.bam} \
            | samtools mpileup -A -B -m {params.indel_min_read} -C {params.mapq_adj} -p -d {params.max_depth} -f {input.ref} - \
            > {output.pileup}

        java -Xmx${{mem_gb}}G -jar $VARSCAN2_BIN mpileup2snp {output.pileup} \
                --min-var-freq {params.min_freq} --p-value {params.pval} --strand-filter {params.strand} --output-vcf 1 \
                > {output.vcf_snp}
        java -Xmx${{mem_gb}}G -jar $VARSCAN2_BIN mpileup2indel {output.pileup} \
                --min-var-freq {params.min_freq} --p-value {params.pval} --strand-filter {params.strand} --output-vcf 1 \
                > {output.vcf_indel}
        '''

