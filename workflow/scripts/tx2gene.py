#!/usr/bin/python3

### Adapted from Danny Bergeron's code

import pandas as pd
import re

gtf_file = snakemake.input.gtf
out_file = snakemake.output.tsv

def transcript2gene(gtf_file):
    # Creating a transcript -> gene dictionary
    tr_dict = dict()
    with open(str(gtf_file)) as f:
        for line in f:
            trans_id = ''
            gene_id = ''

            if 'gene_id' in line:
                gene_id = re.search(
                    r'gene_id "(.*?)"[;,]', line
                ).group(1)
                gene_id = gene_id.split(',')[0]

            if 'transcript_id' in line:
                trans_id = re.search(
                    r'transcript_id "(.*?)"[;]', line
                ).group(1)

            # Insert in tr dict
            if trans_id and gene_id and trans_id not in tr_dict:
                tr_dict[trans_id] = gene_id

            gene_id = ''
            trans_id = ''
    return tr_dict

# Generate the dictionary
_dict = transcript2gene(gtf_file)

with open(out_file, 'w') as f:
    for key, value in _dict.items():
        f.write(str(key) + '\t' + value + '\n')
