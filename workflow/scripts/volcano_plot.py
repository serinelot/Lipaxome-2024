#!/usr/bin/python3

### Adapted from the script written by Danny Bergeron

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np

# Get list of genes after TPM filtering
genes_file = snakemake.input.filtered_genes
genes = pd.read_csv(genes_file, sep='\t', usecols=['gene'])

# Configure comparison, pval_threshold and colors
comp1, comp2 = str(snakemake.wildcards.comp).split('-')

#DE_tool, quantifier = str(snakemake.wildcards.DE_tool), str(snakemake.wildcards.quant)
pval_threshold = snakemake.params.pval_threshold
#colors = {'|log2FC| > 1 & padj < '+str(pval_threshold): 'blue','n.s.': 'grey'}
colors = {'log2FC < -1 & padj < '+str(pval_threshold): 'cornflowerblue','log2FC > 1 & padj < '+str(pval_threshold): 'firebrick','n.s.': 'grey'}

# Load DE df
df = pd.read_csv(snakemake.input.DE_output)

# Rename columns
df.set_axis(['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'], axis=1, inplace=True)

# Drop genes/transcripts with NaN in log2FoldChange, pvalue and/or padj
df = df.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])

# Drop genes that are not in the TPM filtered list
df = df[df.gene.isin(genes.gene)] 

# Filter genes by biotype
# only keep protein coding genes
gtf = snakemake.input.gtf
df_gtf = read_gtf(gtf)
df = df[df.gene.isin(df_gtf.gene_id)]

# Create -log10padj column
df['-log10padj'] = -np.log10(df['padj'])

# Create hue column for significative points (|log2FC| > 1 & padj<0.05)
df.loc[(df['padj'] < pval_threshold) & (df['log2FoldChange'] > 1),
        'sig.'] = 'log2FC > 1 & padj < '+str(pval_threshold)
df.loc[(df['padj'] < pval_threshold) & (df['log2FoldChange'] < -1),
        'sig.'] = 'log2FC < -1 & padj < '+str(pval_threshold)
df['sig.'] = df['sig.'].fillna('n.s.')

# Extract genes that are significant
outfile_up = snakemake.output.up_genes
outfile_down = snakemake.output.down_genes
df_genes = df[~df['sig.'].str.contains('n.s.')]
df_genes = df_genes[['gene','log2FoldChange','padj']]

up = df_genes[df_genes['log2FoldChange']>0]
up.sort_values('log2FoldChange',inplace=True,ascending=False)
down = df_genes[df_genes['log2FoldChange']<0]
down.sort_values('log2FoldChange',inplace=True,ascending=False)
up.to_csv(outfile_up, sep='\t', index=False)
down.to_csv(outfile_down, sep='\t', index=False)

# Create volcano function
def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            pval_threshold, **kwargs):
    """
    Creates a volcano plot (using a x, y and hue column).
    """
    #sns.set_theme()
    plt.figure(figsize=(5.5,4))
    plt.rcParams['svg.fonttype'] = 'none'

    plt.suptitle(title, fontsize=16)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    # Add threshold lines (padj)
    plt.axhline(y=-np.log10(pval_threshold), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(2), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(0.5), color='black', ls='--', lw=0.5)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    leg = plt.legend(fontsize="medium",bbox_to_anchor=(1.6, 0.4), loc='lower right', ncol=1,
        facecolor='white',framealpha=1)
    leg.get_frame().set_linewidth(0)
    plt.savefig(path, bbox_inches='tight', dpi=600)
    
    return

# Create volcano
volcano(df, 'log2FoldChange', '-log10padj', 'sig.',
        'log2(Fold Change)', '-log10(FDR-adjusted p-value)',
        f'{comp1} vs {comp2} using DESeq2',
        colors, snakemake.output.volcano, pval_threshold)
