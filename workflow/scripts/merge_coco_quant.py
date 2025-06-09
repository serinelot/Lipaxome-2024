#!/usr/bin/env python3
"""
Fusionne les sorties CoCo (count, CPM, TPM) de tous les échantillons
en trois matrices consolidées.
Entrées : un TSV par échantillon avec colonnes gene_id, gene_name, count, cpm, tpm
Sorties :
  - merged_counts.tsv : matrice des counts bruts
  - merged_cpm.tsv    : matrice des CPM
  - merged_tpm.tsv    : matrice des TPM
"""
import pandas as pd
from pathlib import Path

# Récupération des chemins depuis Snakemake
in_files   = snakemake.input.coco_tsv
out_counts = Path(snakemake.output.merged_counts)
out_cpm    = Path(snakemake.output.merged_cpm)
out_tpm    = Path(snakemake.output.merged_tpm)

# Création du dossier de sortie si nécessaire
out_dir = out_counts.parent
out_dir.mkdir(parents=True, exist_ok=True)

# Initialisation des DataFrames
counts_df = None
cpm_df    = None
tpm_df    = None

# Parcours de chaque fichier d’entrée
for f in in_files:
    df = pd.read_csv(f, sep="\t")
    sample = Path(f).stem  # ex. "LipC03"

    # Sélection et renommage des colonnes pour ce sample
    df_counts = df[["gene_id","gene_name","count"]].rename(columns={"count": sample})
    df_cpm    = df[["gene_id","gene_name","cpm"]].rename(columns={"cpm": sample})
    df_tpm    = df[["gene_id","gene_name","tpm"]].rename(columns={"tpm": sample})

    # Premier échantillon ? on initialise les tables
    if counts_df is None:
        counts_df = df_counts
        cpm_df    = df_cpm
        tpm_df    = df_tpm
    else:
        # Sinon on merge outer pour garder tous les gènes
        counts_df = counts_df.merge(df_counts, on=["gene_id","gene_name"], how="outer")
        cpm_df    = cpm_df.merge(df_cpm,     on=["gene_id","gene_name"], how="outer")
        tpm_df    = tpm_df.merge(df_tpm,     on=["gene_id","gene_name"], how="outer")

# Écriture des fichiers de sortie
counts_df.to_csv(out_counts, sep="\t", index=False)
cpm_df.   to_csv(out_cpm,    sep="\t", index=False)
tpm_df.   to_csv(out_tpm,    sep="\t", index=False)
