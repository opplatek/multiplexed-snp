rule adapter_scan:
    input:
        "data/raw/{run}/{sample}_L001_R{read}_001.fastq.gz",
    output:
        minion=temp(OUTDIR + "{run}/qc/adapters/{sample}_L001_R{read}_001.minion.fasta"),
        swan=OUTDIR + "{run}/qc/adapters/{sample}_L001_R{read}_001.minion.txt",
    conda:
        "../envs/qc.yaml"
    log:
        OUTDIR + "{run}/qc/adapters/logs/{sample}_L001_R{read}_001.minion.log"
    params:
        adapters_db="db/adapters_merge.txt",
    shell:
        '''
        minion search-adapter -i {input} -show 3 -write-fasta {output.minion} &> {log}
        swan -r {params.adapters_db} -q {output.minion} > {output.swan}
        '''

# TODO: Make match in params
rule adapter_scan_collect:
    input:
        expand(OUTDIR + "{run}/qc/adapters/{sample}_L001_R{read}_001.minion.txt", run=RUNS, sample=SAMPLES, read=READS),
    output:
        txt=OUTDIR + "{run}/qc/adapters/adapters.txt",
        fasta=OUTDIR + "{run}/qc/adapters/adapters.fa"
    shell:
        '''
        match=$(printf '|%.0s' {{1..12}}) # Number of matches between found overrepresented sequence and an adapter in the db

        for i in {input}; do
            echo $(basename $i .minion.txt) >> {output.txt}

            echo -ne ">" >> {output.fasta}
            echo $(basename $i .minion.txt) >> {output.fasta}

            if grep -q "$match" $i; then # If grep doesn't find a match it reports a hidden fail and snakemake rule fails. We have to scan for match before we even attempt the match
                grep -F -B1 -A1 $match $i | sed 's/^ //g' >> {output.txt}
                grep -F -A1 $match $i | sed 's/^ //g' | awk '{{print $1}}' | sed '/^||||||||||||/d' >> {output.fasta}
            fi
        done
        '''
