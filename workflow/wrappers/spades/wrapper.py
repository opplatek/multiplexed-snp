__author__ = "Jan Oppelt"
__copyright__ = "Copyright 2023, Jan Oppelt"
__email__ = "jan.oppelt@gmail.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

input = snakemake.input
output = snakemake.output
threads = snakemake.threads
log = snakemake.log

extra = snakemake.params.get("extra", "")

if hasattr(snakemake.resources, "mem_mb"):
    mem_gb = snakemake.resources.mem_mb // 1000
    memory_requirements = f" --memory {mem_gb}"
else:
    memory_requirements = ""

output_dir = os.path.split(snakemake.output.contigs)[0]

shell("spades.py --careful --threads {threads} {memory_requirements} "
"--tmp-dir {output.dir} "
"-k 21,33,55,77,99,127 -1 {input.r1} -2 {input.r2} -o {output_dir} {extra} &> {log}")
