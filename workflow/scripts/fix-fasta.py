#!/usr/bin/env python3
#
# Fix problematic fasta sequences
#
# * remove spaces in fasta name
# * replace problematic characters
# * split multi-fasta to single-fasta
# * convert multi-line to single-line fasta sequence
# * remove empty lines
# * replace non-IUPAC (ATUGCatugc) characters with N
#

import sys
import os
import re

def split_fasta(filename, odir):
    with open(filename) as f:
        lines = f.readlines()

    header = '' # Init header
    seq = '' # Init seq
    char_replace = f"\s|\.|\||:|-|\+" # Characters that can break the processing to be replaced by "_"
    char_remove = f",|_$" # Characters to completely remove 

    for line in lines:
        if line[0] == '>': # Check for header
            if header:
                with open(f'{odir}/{header}.fa', 'w') as outfile:
                    outfile.write(">" + header + "\n" + seq)

            header = line[1:].strip()
            header = re.sub(char_replace, "_", header)
            header = re.sub(char_remove, "", header)
#            header = re.sub(re.compile(r"(_)\1+"), "_", header) # Replace one or more occurence of the "_" character
#            header = re.sub(re.compile(r"(_)\1{1,}"), "_", header) # Replace one or more occurences of the "_" character in the group \1
            header = re.sub('__+', "_", header) # Replace one or more occurences of the "_" character in the group \1

            seq = '' # Reset seq if new header is encountered
        else:
            if line.strip(): # Skip empty lines
                seq += re.sub(r'[^ATUGCatugc]', "N", line.strip()) # Replace non-IUPAC characters with N and add to seq
    
    with open(f'{odir}/{header}.fa', 'w') as outfile: # Write the last sequence
        outfile.write(f">{header}\n{seq}\n")


if len(sys.argv) != 3:
    print("Exactly two arguments (input fasta file and output directory) has to be speficied.")
    sys.exit()
else:
    ifile = sys.argv[1]
    odir = sys.argv[2]
    os.makedirs(odir, exist_ok = True)

    split_fasta(ifile, odir)
