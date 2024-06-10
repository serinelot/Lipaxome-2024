#!/usr/bin/python3

import pandas as pd
from gtfparse import read_gtf
import os

tx2gene = snakemake.input.tx2gene
gtf = snakemake.input.gtf
outfile = snakemake.output.tpm
log_file_path = snakemake.log[0]

final_df = pd.DataFrame()
cycle = 1
sample_list = []

with open(log_file_path, 'w') as log_file:
    log_file.write("Starting merge process\n")

    log_file.write(f"Reading GTF file from {gtf}\n")
    df_gtf = read_gtf(gtf)
    protein_coding_genes = set(df_gtf[df_gtf['gene_biotype'] == 'protein_coding']['gene_id'].dropna().unique())
    protein_coding_transcripts = set(df_gtf[df_gtf['gene_id'].isin(protein_coding_genes)]['transcript_id'].dropna().unique())
    log_file.write(f"Number of protein coding transcripts (based on gene_biotype): {len(protein_coding_transcripts)}\n")

    for q in snakemake.input.quant:
        log_file.write(f"Processing file: {q}\n")
        
        data = pd.read_csv(q, sep='\t', usecols=['target_id', 'tpm'])
        log_file.write(f"Read {data.shape[0]} rows from {q}\n")

        data = data[data['target_id'].isin(protein_coding_transcripts)]
        log_file.write(f"Filtered {data.shape[0]} rows to protein-coding transcripts\n")

        sample = os.path.basename(os.path.dirname(q))

        data.set_index('target_id', inplace=True)
        data.rename(columns={"tpm": sample}, inplace=True)

        if cycle == 1:
            final_df = data
        else:
            final_df = pd.merge(final_df, data, left_index=True, right_index=True, how='outer')

        sample_list += [sample]
        cycle += 1

        log_file.write(f"Final dataframe shape after merging {sample}: {final_df.shape}\n")

    log_file.write(f"Reading transcript to gene ID mapping from {tx2gene}\n")
    ids = pd.read_csv(tx2gene, sep='\t', names=['transcript', 'gene'])
    ids.set_index('transcript', inplace=True, drop=False)

    final_df = pd.merge(final_df, ids, left_index=True, right_index=True, how='left')
    final_df.set_index('gene', inplace=True)
    log_file.write(f"Final dataframe shape after adding gene IDs: {final_df.shape}\n")

    # Skipping TPM filtering to keep all rows
    log_file.write("Skipping TPM filtering\n")
    filtered = final_df
    log_file.write(f"Dataframe shape without TPM filtering: {filtered.shape}\n")

    id_name = df_gtf[['gene_id', 'gene_name']].drop_duplicates(ignore_index=True)

    index = filtered.index.tolist()
    names = []

    for i in range(len(index)):
        id = index[i]
        if id in id_name['gene_id'].values:
            name = id_name[id_name['gene_id'] == id].iloc[0]['gene_name']
        else:
            name = "Unknown"
        names.append(name)

    log_file.write("Adding gene names\n")
    filtered['gene_name'] = names
    cols = filtered.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    filtered = filtered[cols]

    log_file.write("Writing final output file\n")
    filtered.to_csv(outfile, sep='\t')
    log_file.write("Merge process completed successfully\n")
