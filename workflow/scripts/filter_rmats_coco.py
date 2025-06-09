#!/usr/bin/env python3
import pandas as pd
import os
import shutil
import warnings
import sys
from pathlib import Path

# — Snakemake inputs & outputs —
summary_file = snakemake.input.summary
tpm_matrix   = snakemake.input.tpm_matrix
gtf_corr     = snakemake.input.gtf_corr
out_file     = snakemake.output.filtered
log_path     = Path(snakemake.log[0])

fdr     = snakemake.params.fdr
dpsi    = snakemake.params.dpsi
min_tpm = snakemake.params.min_tpm

def read_gtf(filepath):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df = pd.read_csv(
            filepath,
            sep="\t",
            comment="#",
            header=None,
            on_bad_lines="warn"
        )
        df.columns = [
            "seqname","source","feature","start","end",
            "score","strand","frame","attributes"
        ]
        for warn in w:
            print(f"Warning: {warn.message}", file=sys.stderr)
        return df

# Création du répertoire filtered
raw_dir = os.path.dirname(summary_file)
filtered_dir = os.path.dirname(out_file)
os.makedirs(filtered_dir, exist_ok=True)

with log_path.open("w") as log:
    log.write("=== Filtrage rMATS (CoCo‑corrigée) ===\n")
    log.write(f"Chargement TPM matrix: {tpm_matrix}\n")
    tpm_df = pd.read_csv(tpm_matrix, sep="\t", index_col=0)

    log.write(f"Lecture du GTF corrigé: {gtf_corr}\n")
    df_gtf = read_gtf(gtf_corr)
    # extraire (gene_id,gene_biotype)
    gid = df_gtf["attributes"].str.extract(r'gene_id "([^"]+)"')[0]
    biotype = df_gtf["attributes"].str.extract(r'gene_biotype "([^"]+)"')[0]
    df_biotype = pd.DataFrame({"gene_id": gid, "gene_biotype": biotype})
    pc_genes = set(df_biotype.query("gene_biotype=='protein_coding'")["gene_id"].unique())
    log.write(f"→ {len(pc_genes)} gènes protein_coding détectés dans GTF CoCo\n")

    def filter_threshold(df):
        return df[(df["FDR"] <= fdr) & (df["IncLevelDifference"].abs() >= dpsi)]

    def filter_tpm(df):
        # on ne garde que les gènes exprimés au-dessus du min_tpm en TPM
        expressed = set(tpm_df.index[tpm_df.max(axis=1) >= min_tpm])
        return df[df["GeneID"].isin(expressed & pc_genes)]

    events = [
        "SE.MATS.JC.txt",
        "A5SS.MATS.JC.txt",
        "A3SS.MATS.JC.txt",
        "MXE.MATS.JC.txt",
        "RI.MATS.JC.txt"
    ]

    for ev in events:
        inpath = os.path.join(raw_dir, ev)
        if not os.path.exists(inpath):
            print(f"⚠️  Manque {inpath}", file=sys.stderr)
            continue
        log.write(f"Processing {ev}\n")
        df = pd.read_csv(inpath, sep="\t")
        df = filter_threshold(df)
        df = filter_tpm(df)
        prefix = ev.split(".")[0]
        outpath = os.path.join(filtered_dir, f"{prefix}.tsv")
        df.to_csv(outpath, sep="\t", index=False)
        log.write(f"→ écrit {outpath}\n")

    log.write("=== Filtrage CoCo terminé ===\n")

# Cleanup du tmp si besoin
tmp_dir = os.path.join(raw_dir, "..", "tmp")
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)
