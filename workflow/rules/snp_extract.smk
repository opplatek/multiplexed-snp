rule extract_gatk:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.phased.vcf",
    output:
        temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.gatk.tmp"),
    shell:
        '''
		cat {input} | awk '$7 == "PASS" {{print $1, $2}}' > {output}
        '''


rule extract_samtools:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/samtools/{sample}.recal.samtools.filt.vcf",
    output:
        temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.samtools.tmp"),
    shell:
        '''
        cat {input} | awk '{{print $1, $2}}' | awk '/POS/{{i++}}i' | sed '1d' > {output}
        '''


rule extract_varscan2:
    input:
        vcf_snp=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.snp.vcf",
        vcf_indel=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.indel.vcf",
    output:
        temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.varscan2.tmp"),
    threads: 4
    shell:
        '''
		sort -k2,2n -T {resources.tmpdir} --parallel={threads} \
            <(cat {input.vcf_snp} | awk '{{print $1, $2}}' | grep -v "^#" | sed '1d') \
            <(cat {input.vcf_indel} | awk '{{print $1, $2, "\\tindel"}}' | grep -v "^#" | sed '1d') > {output}
        '''


rule extract_merge:
    input:
        gatk=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.gatk.tmp",
        samtools=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.samtools.tmp",
        varscan2=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.varscan2.tmp",
    output:
        temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/extract/{sample}.extract.snp-raw"),
    shell:
        '''
		echo -e "gatk\\tsamtools\\tvarscan\\tvarscanIndel" > {output}
		paste -d'\\t' {input.gatk} {input.samtools} {input.varscan2} >> {output}
        '''


rule vcf_merge:
    input:
        gatk=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gatk/{sample}.recal.gatk.filt.phased.vcf",
        samtools=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/samtools/{sample}.recal.samtools.filt.vcf",
        varscan2_snp=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.snp.vcf",
        varscan2_indel=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/varscan2/{sample}.recal.varscan2.indel.vcf",
    output:
        gatk_gz=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.gatk.filt.phased.vcf.gz"),
        gatk_gz_tbi=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.gatk.filt.phased.vcf.gz.tbi"),
        samtools_gz=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.samtools.filt.vcf.gz"),
        samtools_gz_tbi=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.samtools.filt.vcf.gz.tbi"),
        varscan2_snp_gz=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.varscan2.snp.vcf.gz"),
        varscan2_snp_gz_tbi=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.varscan2.snp.vcf.gz.tbi"),
        varscan2_indel_gz=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.varscan2.indel.vcf.gz"),
        varscan2_indel_gz_tbi=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.recal.varscan2.indel.vcf.gz.tbi"),
        merged=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/vcf/{sample}.gatkSamtoolsVarScan2_Merged.vcf",
    conda:
        "../envs/snp.yaml"
    shell:
        '''
		sed -i "s/FORMAT\\t{wildcards.sample}/FORMAT\\t{wildcards.sample}_gatk/" {input.gatk}
		sed -i "s/FORMAT\\t{wildcards.sample}/FORMAT\\t{wildcards.sample}_samtools/" {input.samtools}
		sed -i "s/FORMAT\\tSample1/FORMAT\\t{wildcards.sample}_varscanSnp/" {input.varscan2_snp}
		sed -i "s/FORMAT\\tSample1/FORMAT\\t{wildcards.sample}_varscanIndel/" {input.varscan2_indel}

		bgzip -c {input.gatk} > {output.gatk_gz} && tabix {output.gatk_gz}
		bgzip -c {input.samtools} > {output.samtools_gz} && tabix {output.samtools_gz}
		bgzip -c {input.varscan2_snp} > {output.varscan2_snp_gz} && tabix {output.varscan2_snp_gz}
		bgzip -c {input.varscan2_indel} > {output.varscan2_indel_gz} && tabix {output.varscan2_indel_gz}

		vcf-merge \
            {output.gatk_gz} \
            {output.samtools_gz} \
            {output.varscan2_snp_gz} \
            {output.varscan2_indel_gz} > {output.merged}
        '''

rule extract_genotypes:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/vcf/{sample}.gatkSamtoolsVarScan2_Merged.vcf",
    output:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/gt/{sample}.gatkSamtoolsVarScan2_Merged.gt.txt",
    conda:
        "../envs/snp.yaml"
    shell:
        "vcf-to-tab < {input} > {output}"


rule extract_overlaps:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/extract/{sample}.extract.snp-raw",
    output:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/extract/{sample}.extract.snp",
    conda:
        "../envs/snp.yaml"
    script:
        "../scripts/intersect_nonpooled.R"


rule extract:
    input:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/extract/{sample}.extract.snp",
    output:
        OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/snp/{sample}.all.snp",
    shell:
       """
        cat {input} \
        | sed 's/ /\\t/g' \
		| sed 's/\\tindel/indel/g' \
		| sed '1d' \
		| sed '1i nrow\\tgatk\\tgatk\\tsamtools\\tsamtools\\tvarscan\\tvarscan\\tvarscanIndel\\toverlap_all\\toverlap_all\\toverlap_two\\toverlap_two\\n' \
        | awk '{{print $11, $12}}' \
        | tail -n+2 \
        | sort -nk 2,2 -k 1,1 | uniq \
        | sed 's/- -//g' \
        | sed '/^\s*$/d' \
        | sort -k1,1 -k2n \
        > {output}
       """

