#!/usr/bin/python3
import pandas as pd

""" Generate a merged dataframe of all samples output from coco cc (one df for
    the count, one for the cpm and one for the tpm) """

input_counts_paths = snakemake.input.counts
output_counts = snakemake.output.merged_counts
output_cpm = snakemake.output.merged_cpm
output_tpm = snakemake.output.merged_tpm

# Generate the count, cpm and tpm list of sample dfs
count_temp, cpm_temp, tpm_temp = [], [], []
sorted_paths = sorted(input_counts_paths)
for i, file in enumerate(sorted_paths):
    file_name = file.split('/')[-1].split('.')[0]
    df = pd.read_csv(file, sep='\t')
    if i == 0:  # Add gene_id and gene_name column in the list only one time (for the first sample here)
        gene_id = df[['gene_id', 'gene_name']]
        count_temp.append(gene_id)
        cpm_temp.append(gene_id)
        tpm_temp.append(gene_id)

    count, cpm, tpm = df[['count']], df[['cpm']], df[['tpm']]
    count.columns, cpm.columns, tpm.columns = [file_name], [file_name], [file_name]
    count_temp.append(count)
    cpm_temp.append(cpm)
    tpm_temp.append(tpm)

# Generate the count, cpm and tpm dataframes
count_df, cpm_df, tpm_df = pd.concat(count_temp, axis=1), pd.concat(cpm_temp, axis=1), pd.concat(tpm_temp, axis=1)
count_df.to_csv(output_counts, index=False, sep='\t')
cpm_df.to_csv(output_cpm, index=False, sep='\t')
tpm_df.to_csv(output_tpm, index=False, sep='\t')
