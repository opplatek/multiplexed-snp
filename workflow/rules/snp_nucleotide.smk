rule nucleotide_snp:
    input:
        bam=OUTDIR + "{run}/" + OUTDIR_SUB + "bams/{sample}.duplRem.bam",
        ref="data/references/{run}/{sample}_ref.fa",
        snp=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/snp/{sample}.all.snp",
    output:
        tmp=temp(OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/{sample}.snp.bamreadcount.tmp"),
        full=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/snp/{sample}.full.nucDist.txt",
        long=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/snp/{sample}.long.nucDist.txt",
        short=OUTDIR + "{run}/" + OUTDIR_SUB + "variants/{sample}/snp/{sample}.short.nucDist.txt",
    conda:
        "../envs/snp.yaml"
    params:
        insert_line=r"chr	position	reference_base	depth	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end",
        remove_line=r"\t=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00",
        min_mapq=10,
        min_phred=15,
    shell:
        '''
        awk '{{print $1, $2, $2}}' {input.snp} > {output.tmp}

        bam-readcount -q {params.min_mapq} -b {params.min_phred} -l {output.tmp} -w 1 -f {input.ref} {input.bam} \
        | sed "1s/^/{params.insert_line}\\n/" \
        | sed "s/{params.remove_line}//g" \
        | sed 's/:/\\t/g' > {output.full}

        awk 'BEGIN {{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $10, $11, $19, $20, $24, $25, $33, $34, $38, $39, $47, $48, $52, $53}}' {output.full} > {output.long}
        awk 'BEGIN {{OFS="\\t"}}{{print $1, $2, $3, $4, $5, $6, $19, $20, $33, $34, $47, $48}}' {output.full} > {output.short}
        '''


