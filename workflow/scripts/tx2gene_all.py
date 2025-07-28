#!/usr/bin/env python3
import re
from pathlib import Path

# Snakemake injecte ici :
#   snakemake.input.gtf
#   snakemake.output.tx2gene_all
gtf_path   = Path(snakemake.input.gtf)
out_path   = Path(snakemake.output.tx2gene_all)

tx2gene = {}

with gtf_path.open("rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols    = line.rstrip("\n").split("\t")
        feature = cols[2]
        # Certains GTF n'ont pas de "transcript", on prend aussi les exon
        if feature not in ("transcript", "exon"):
            continue

        attrs = cols[8]
        m_gid = re.search(r'gene_id "([^"]+)"', attrs)
        m_tid = re.search(r'transcript_id "([^"]+)"', attrs)
        if not m_gid or not m_tid:
            continue
        gid = m_gid.group(1)
        tid = m_tid.group(1)
        # ne pas écraser si déjà rempli
        if tid not in tx2gene:
            tx2gene[tid] = gid

# Création du dossier si besoin
out_path.parent.mkdir(parents=True, exist_ok=True)

# Écriture du TSV avec en-tête
with out_path.open("w") as out:
    out.write("transcript_id\tgene_id\n")
    for tid in sorted(tx2gene):
        out.write(f"{tid}\t{tx2gene[tid]}\n")
