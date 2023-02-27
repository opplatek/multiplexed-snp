# TODO:
# rule ref_split:
#     input:

#     output:
#         dir()
#     shell:
#         "python3 ../scripts/fix-fasta.py {input} {output}"


# rule ref_merge:
#     input:
#         from panda df
#     output:
#         "data/references/{run}/{sample}_ref.fa",
#     shell:
#         "cat {input} > {output}"


# rule ref_fix:
#     input:
#         from panda df
#     output:
#         "data/references/{run}/{sample}_ref.fa",
#     shell:
#         "python3 ../scripts/fix-fasta.py {input} {output}"


rule index_dos2unix:
    input:
        "data/references/{run}/{sample}_ref.fa",
    output:
        temp(touch("data/references/{run}/{sample}_ref.fa.prep")),
    shell:
        '''
        dos2unix {input} # Make sure the reference has UNIX formating
        '''

