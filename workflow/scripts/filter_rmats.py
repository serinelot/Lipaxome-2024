#!/usr/bin/env python3

import pandas as pd
import os
import shutil
import warnings
import sys

# Chemin vers le répertoire temporaire où snakemake stocke les fichiers intermédiaires
tmp_dir = snakemake.params.dir + "/tmp"

# Fonction pour lire le fichier GTF en affichant les lignes mal formées
def read_gtf(filepath):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, on_bad_lines='warn')
        df.columns = [
            'seqname', 'source', 'feature', 'start', 'end', 'score',
            'strand', 'frame', 'attributes'
        ]

        # Imprimer les avertissements pour les lignes mal formées
        for warning in w:
            print(f"Warning: {warning.message}", file=sys.stderr)
        
        return df

# Paramètres depuis snakemake
SPLICING_EVENTS = ['SE.MATS.JC.txt', 'A5SS.MATS.JC.txt', 'A3SS.MATS.JC.txt', 'MXE.MATS.JC.txt', 'RI.MATS.JC.txt']
raw_dir = snakemake.params.dir + "/raw"
out_dir = snakemake.params.dir + "/filtered"
tpm = snakemake.params.tpm
fdr = snakemake.params.fdr
deltapsi = snakemake.params.deltapsi
gtf = snakemake.params.gtf

# Lecture de la matrice TPM
tpm_df = pd.read_csv(tpm, sep='\t')

# Lecture du fichier GTF
df_gtf = read_gtf(gtf)
# Filtrage pour obtenir les gènes codant pour des protéines
id_biotype = df_gtf['attributes'].str.extract('gene_id "([^"]+)"')[0]
id_biotype = pd.concat([id_biotype, df_gtf['attributes'].str.extract('gene_biotype "([^"]+)"')[0]], axis=1)
id_biotype.columns = ['gene_id', 'gene_biotype']
pc_genes_list = id_biotype[id_biotype['gene_biotype'] == 'protein_coding']['gene_id'].unique().tolist()

def filter_by_threshold(df, fdr, deltapsi):
    """Filtre les événements d'épissage en fonction du FDR et de la différence de niveau d'inclusion."""
    return df[(df['FDR'] <= fdr) & (df['IncLevelDifference'].abs() >= deltapsi)]

def filter_by_tpm(rmats_df, tpm_df, pc_genes_list):
    """Filtre pour conserver les gènes exprimés dans au moins un échantillon."""
    exp_genes = tpm_df[tpm_df.max(axis=1) >= 1]['gene'].unique()
    pc_genes_set = set(pc_genes_list)
    return rmats_df[rmats_df['GeneID'].isin(exp_genes) & rmats_df['GeneID'].isin(pc_genes_set)]

# Traitement pour chaque type d'événement d'épissage
for event in SPLICING_EVENTS:
    file_path = os.path.join(raw_dir, event)
    df = pd.read_csv(file_path, sep='\t')
    df = filter_by_threshold(df, fdr, deltapsi)
    df = filter_by_tpm(df, tpm_df, pc_genes_list)
    # Sauvegarde des événements d'épissage filtrés
    event_type = event.split('.')[0]
    df.to_csv(os.path.join(out_dir, event_type + '.tsv'), sep='\t', index=False)

# Nettoyage du répertoire temporaire
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)