# TODO: Add automatic detection of FASTA/FASTQ
# TODO: Add automatic detection of PE/SE reads

# Decide if we have global variable MIN_LEN_FASTQ for fastq length filtering and choose the output - w/ or w/o length filtering
def decide_fastq_input_lenfilt(wildcards):
    if 'MIN_LEN_FASTQ' in globals() and MIN_LEN_FASTQ > 1:
        return {'r1': "data/prep/{run}/{sample}_L001_R1_001_lenfilt.fastq.gz",'r2': "data/prep/{run}/{sample}_L001_R2_001_lenfilt.fastq.gz"}
    else:
        return {'r1': "data/prep/{run}/{sample}_L001_R1_001_trim.fastq.gz",'r2': "data/prep/{run}/{sample}_L001_R2_001_trim.fastq.gz"}

# TODO: Move this into a Python wrapper for easier tracking of what command has been used
rule trim_ends:
    input:
        "data/raw/{run}/{sample}_L001_R{read}_001.fastq.gz",
    output:
        temp("data/prep/{run}/{sample}_L001_R{read}_001_trimends.fastq.gz"),
    conda:
        "../envs/preprocess.yaml"         
    params:
        trim5=5,
        trim3=0,
    shell:
        '''
        if [[ {params.trim5} > 0 && {params.trim3} > 0 ]]; then
            echo "Trimming both ends of {input} (seqtk)."
            seqtk trimfq -b {params.trim5} -e {params.trim3} {input} | gzip -c > {output}
        elif [ {params.trim5} > 0 ]; then
            echo "Trimming 5' end of {input} (seqtk)."
            seqtk trimfq -b {params.trim5} {input} | gzip -c > {output}
        elif [ {params.trim3} > 0 ]; then
            echo "Trimming 3' end of {input} (seqtk)."
            seqtk trimfq -e {params.trim3} {input} | gzip -c > {output}
        else
            echo "Not trimming anything of {input} (seqtk)."
            ln -sf {input} {output}
        fi
        '''


rule preprocess_trimmomatic_pe:
    input:
        r1="data/prep/{run}/{sample}_L001_R1_001_trimends.fastq.gz",
        r2="data/prep/{run}/{sample}_L001_R2_001_trimends.fastq.gz",
    output:
        r1_pair="data/prep/{run}/{sample}_L001_R1_001_trim.fastq.gz",
        r2_pair="data/prep/{run}/{sample}_L001_R2_001_trim.fastq.gz",
        r1_unpair=temp("data/prep/{run}/{sample}_L001_R1_001_trim_unpaired.fastq.gz"),
        r2_unpair=temp("data/prep/{run}/{sample}_L001_R2_001_trim_unpaired.fastq.gz"),
    conda:
        "../envs/preprocess.yaml"        
    log:
        OUTDIR + "{run}/qc/preprocess/logs/{sample}.trimmomatic.log",
    threads: 32
    params:
        phred=20,
        read_len_trim=250,
        min_len=35,
    shell:
        '''
        trimmomatic PE -phred33 -threads {threads} {input.r1} {input.r2} \
            {output.r1_pair} {output.r1_unpair} \
            {output.r2_pair} {output.r2_unpair} \
            CROP:{params.read_len_trim} \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.phred} \
            MINLEN:{params.min_len} &> {log}
        '''


rule filter_length:
    input:
        r1="data/prep/{run}/{sample}_L001_R1_001_trim.fastq.gz",
        r2="data/prep/{run}/{sample}_L001_R2_001_trim.fastq.gz",
    output:
        r1_filt="data/prep/{run}/{sample}_L001_R1_001_lenfilt.fastq.gz",
        r2_filt="data/prep/{run}/{sample}_L001_R2_001_lenfilt.fastq.gz",
    conda:
        "../envs/preprocess.yaml"        
    params:
        mem_mb=2000,
        min_len=MIN_LEN_FASTQ if 'MIN_LEN_FASTQ' in globals() else 1,
    shell:
        '''
        reformat.sh -Xmx{params.mem_mb}m in={input.r1} in2={input.r2} \
            out={output.r1_filt} out2={output.r2_filt} \
            minlength={params.min_len}
        '''


rule link_preprocess_pe:
    input:
        unpack(decide_fastq_input_lenfilt)
    output:
        r1_prep="data/prep/{run}/{sample}_L001_R1_001_prep.fastq.gz",
        r2_prep="data/prep/{run}/{sample}_L001_R2_001_prep.fastq.gz",
    shell:
        '''
        ln -sf `basename {input.r1}` {output.r1_prep}
        ln -sf `basename {input.r2}` {output.r2_prep}
        '''

