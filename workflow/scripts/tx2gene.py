#!/usr/bin/env python3
"""
Génère un fichier tx2gene.tsv à partir d'un fichier GTF.

Pour chaque transcript_id du GTF, on extrait la gene_id associée.
Le résultat est écrit en deux colonnes : transcript_id \t gene_id
"""

import re
from pathlib import Path

# Récupération des chemins via Snakemake
gtf_path = Path(snakemake.input.gtf)
out_path = Path(snakemake.output.tsv)

# Dictionnaire transcript → gene
tx2gene = {}

with gtf_path.open("r") as f:
    for line in f:
        # Ignorer les commentaires
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        feature = cols[2]

        # On ne garde que les transcripts (ou exons si besoin)
        if feature != "transcript":
            continue

        attrs = cols[8]

        # Recherche des IDs
        gene_m = re.search(r'gene_id "([^"]+)"', attrs)
        tx_m   = re.search(r'transcript_id "([^"]+)"', attrs)

        if not gene_m or not tx_m:
            continue

        gene_id      = gene_m.group(1)
        transcript_id = tx_m.group(1)

        # Ne pas écraser si déjà trouvé (première occurrence)
        if transcript_id not in tx2gene:
            tx2gene[transcript_id] = gene_id

# Écriture du fichier de sortie, trié par transcript_id
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as f:
    for tx in sorted(tx2gene):
        f.write(f"{tx}\t{tx2gene[tx]}\n")
