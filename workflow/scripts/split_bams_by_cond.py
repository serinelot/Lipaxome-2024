#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# 1) Lire le design
design = pd.read_csv("data/design.tsv", sep="\t")

# 2) Chemin source des BAM primaires
src_base = Path("results/splicing/star")

# 3) Répertoires de destination par condition
dest_base = Path("results/splicing/star")
dest_dirs = {
    "FXS":     dest_base / "FXS",
    "Control": dest_base / "Control",
}

# 4) Ouvrir le script à générer
script_path = Path("scripts/split_bams_by_cond.sh")
script_path.parent.mkdir(parents=True, exist_ok=True)

with script_path.open("w") as out:
    out.write("#!/usr/bin/env bash\n\n")
    # Créer les répertoires de destination
    for d in dest_dirs.values():
        out.write(f"mkdir -p {d}\n")
    out.write("\n")
    # Pour chaque échantillon, copier le BAM & son index .bai au bon endroit
    for _, row in design.iterrows():
        samp      = row["sample"]
        cond      = row["condition"]
        src_dir   = src_base / samp
        bam_file  = src_dir / f"{samp}_Aligned.sortedByCoord.out.primary.bam"
        bai_file  = src_dir / f"{samp}_Aligned.sortedByCoord.out.primary.bam.bai"
        dest_dir  = dest_dirs[cond]
        dest_bam  = dest_dir / bam_file.name
        dest_bai  = dest_dir / bai_file.name

        out.write(f"cp {bam_file} {dest_bam}\n")
        out.write(f"cp {bai_file} {dest_bai}\n")
