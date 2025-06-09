#!/usr/bin/env python3
"""
Fusionne les abundance.tsv de Kallisto en une matrice TPM
pour les gènes protéinocodants.
"""

import pandas as pd
from gtfparse import read_gtf
from pathlib import Path

# Snakemake inputs
quant_files  = snakemake.input.abundances
gtf_path     = snakemake.input.gtf
tx2gene_tsv  = snakemake.input.tx2gene
out_path     = Path(snakemake.output.tpm_matrix)
log_path     = Path(snakemake.log[0])

# Préparer le dossier de sortie
out_path.parent.mkdir(parents=True, exist_ok=True)

with log_path.open("w") as log:
    log.write("=== Merge kallisto TPM quantifications ===\n")
    log.write(f"Reading filtered GTF: {gtf_path}\n")
    df_gtf = read_gtf(gtf_path)

    # Extraire les transcripts protéinocodants
    pc_genes = (
        df_gtf.loc[df_gtf["gene_biotype"] == "protein_coding", "gene_id"]
        .dropna()
        .unique()
    )
    pc_transcripts = (
        df_gtf.loc[df_gtf["gene_id"].isin(pc_genes), "transcript_id"]
        .dropna()
        .unique()
    )
    log.write(f"Found {len(pc_transcripts)} protein-coding transcripts\n")

    merged = None
    for q in quant_files:
        log.write(f"  * Loading {q}\n")
        # On ne prend que la colonne TPM
        df = pd.read_csv(q, sep="\t", usecols=["target_id", "tpm"])
        df = df[df["target_id"].isin(pc_transcripts)]
        sample = Path(q).parent.name
        log.write(f"    → kept {df.shape[0]} rows for sample {sample}\n")

        df = df.set_index("target_id").rename(columns={"tpm": sample})
        merged = df if merged is None else merged.join(df, how="outer")

    log.write("Merging with tx2gene mapping\n")
    tx2g = (
        pd.read_csv(tx2gene_tsv, sep="\t", names=["transcript","gene"])
        .set_index("transcript")
    )
    merged = merged.join(tx2g, how="left")

    # Reset index and rename the transcript column correctly
    merged = (
        merged
        .reset_index()                        # brings target_id back as a column
        .rename(columns={"target_id": "transcript"})
        .set_index("gene")                   # now switch to gene as index
        .drop(columns=["transcript"])        # drop the transcript column
    )

    log.write("Adding gene names\n")
    names = (
        df_gtf.loc[df_gtf["gene_biotype"] == "protein_coding", ["gene_id","gene_name"]]
        .drop_duplicates()
        .set_index("gene_id")["gene_name"]
    )
    merged["gene_name"] = merged.index.map(lambda g: names.get(g, "NA"))

    # Mettre gene_name en première colonne
    cols = ["gene_name"] + [c for c in merged.columns if c != "gene_name"]
    merged = merged[cols]

    log.write(f"Writing merged TPM matrix to {out_path}\n")
    merged.to_csv(out_path, sep="\t")
    log.write("=== Merge completed successfully ===\n")
