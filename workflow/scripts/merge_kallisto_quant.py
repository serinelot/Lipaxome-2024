#!/usr/bin/python3

import pandas as pd
from gtfparse import read_gtf
import os

tx2gene = snakemake.input.tx2gene
gtf = snakemake.input.gtf
outfile = snakemake.output.tpm

final_df = pd.DataFrame()
cycle = 1
sample_list = []

for q in snakemake.input.quant:

    data = pd.read_csv(q, sep='\t', usecols=['target_id','tpm'])

    # get sample name
    sample = os.path.basename(os.path.dirname(q))

    # reformat dataframe
    data.set_index('target_id', inplace=True)
    data.rename(columns={"tpm":sample}, inplace=True)

    # Merge dataframes
    if cycle == 1:
        final_df = data
    else:
        final_df = pd.merge(final_df, data, left_index=True, right_index=True)
    
    sample_list += [sample]
    cycle+=1

# transcript ID --> gene ID
ids = pd.read_csv(tx2gene, sep='\t', names=['transcript','gene'])
ids.set_index('transcript', inplace=True, drop=False)

final_df = pd.merge(final_df, ids, left_index=True, right_index=True)
final_df.set_index('gene', inplace=True)

# Filter genes by TPM
# TPM >=1 in at least one sample
filtered = final_df[~((final_df[sample_list]<1).all(axis=1))]

# Filter genes by biotype
# only keep protein coding genes
df_gtf = read_gtf(gtf)
filtered_pc = filtered[filtered.index.isin(df_gtf.gene_id)]

# Add gene name
id_name = df_gtf[['gene_id','gene_name']].drop_duplicates(ignore_index=True)

index = filtered_pc.index.tolist()
names =[]

for i in range(len(index)):
    id = index[i]
    name = id_name[id_name['gene_id']==id].iloc[0]['gene_name']
    names.append(name)

print(names)
filtered_pc['gene_name'] = names
cols = filtered_pc.columns.tolist()
cols = cols[-1:] + cols[:-1]
cols = cols[-1:] + cols[:-1]
filtered_pc = filtered_pc[cols]

# Write to file
filtered_pc.to_csv(outfile, sep='\t')
