import pandas as pd
from pathlib import Path

# Chemins depuis config/config.json
fastq_dir   = config["path"]["fastq"].rstrip("/") + "/"       # e.g. "data/fastq/"
sra_list    = config["fastq"]["sra_list"]                     # e.g. "data/references/sra_id_wget.txt"
rename_csv  = config["fastq"]["rename_tsv"]                   # e.g. "data/references/fastq_rename.csv"

# ─── Règle 1 : téléchargement ────────────────────────────────────────────────────
rule fastq_download:
    """
    Télécharger tous les FASTQ listés dans sra_list vers data/fastq/,
    puis marquer la fin avec data/references/download_fastq.done
    """
    input:
        sra = sra_list
    output:
        touch("data/references/download_fastq.done")
    run:
        # Créer le dossier de destination
        Path(fastq_dir).mkdir(parents=True, exist_ok=True)

        # Lancer wget
        shell(f"wget -i {input.sra} -P {fastq_dir}")

        # Vérifier qu’au moins un .fastq.gz est présent
        fastqs = list(Path(fastq_dir).glob("*.fastq.gz"))
        if not fastqs:
            raise RuntimeError(f"Aucun FASTQ téléchargé dans {fastq_dir}")

        # Créer le marqueur
        Path(output[0]).touch()


# ─── Règle 2 : renommage ─────────────────────────────────────────────────────────
rule rename_fastq:
    """
    Renommer les FASTQ selon fastq_rename_fastq.csv (sep=';').
    On saute les anciens fichiers manquants avec un warning.
    """
    input:
        csv = rename_csv,
        done = "data/references/download_fastq.done"   # on attend le download d’abord
    output:
        touch("data/references/rename_fastq.done")
    run:
        from csv import DictReader

        df = pd.read_csv(input.csv, sep=";", dtype=str)
        df = df.dropna(subset=["old_fastq_name","new_fastq_name"])
        df = df[df["old_fastq_name"].str.strip() != ""]
        df = df[df["new_fastq_name"].str.strip() != ""]

        for row in df.itertuples(index=False):
            old = row.old_fastq_name.strip()
            new = row.new_fastq_name.strip()
            src = Path(fastq_dir) / old
            dst = Path(fastq_dir) / new

            if src.exists():
                src.rename(dst)
            else:
                print(f"WARNING: fichier introuvable, saut de {src}")

        # Optionnel : vérifier que **au moins un** renommage a réussi
        renamed = [Path(fastq_dir)/n for n in df.new_fastq_name.str.strip()]
        if not any(p.exists() for p in renamed):
            raise RuntimeError("Aucun fichier n’a été renommé !")

        # Marquer la fin
        Path(output[0]).touch()
