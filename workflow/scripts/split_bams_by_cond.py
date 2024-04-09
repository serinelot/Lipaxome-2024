import pandas as pd

# Lire le fichier de design, en assumant que le répertoire courant est workflow
design = pd.read_csv("data/design.tsv", sep="\t")

# Chemin où se trouvent actuellement les fichiers BAM, relativement à workflow
src_dir = "results/splicing/star/"

# Dossiers de destination basés sur les conditions, relativement à workflow
dest_dirs = {
    "FXS": "results/splicing/star_FXS/",
    "Control": "results/splicing/star_Control/"
}

# Générer les commandes pour créer les dossiers de destination s'ils n'existent pas
with open("scripts/split_bams_files.sh", "w") as f:
    # S'assurer que les dossiers de destination existent
    for dest in dest_dirs.values():
        f.write(f"mkdir -p {dest}\n")
    
    # Générer les commandes de copie des fichiers BAM et leurs fichiers d'index .bai
    for index, row in design.iterrows():
        # Ajuster les chemins source pour correspondre à la structure des fichiers et dossiers
        src_file_bam = f"{src_dir}{row['sample']}/{row['sample']}_Aligned.sortedByCoord.out.primary.bam"
        src_file_bai = f"{src_dir}{row['sample']}/{row['sample']}_Aligned.sortedByCoord.out.primary.bam.bai"
        
        # Définir les chemins de destination, incluant l'identifiant dans le nom du fichier
        dest_file_bam = f"{dest_dirs[row['condition']]}{row['sample']}_Aligned.sortedByCoord.out.primary.bam"
        dest_file_bai = f"{dest_dirs[row['condition']]}{row['sample']}_Aligned.sortedByCoord.out.primary.bam.bai"
        
        # Écrire les commandes pour copier les fichiers BAM et leurs index
        f.write(f"cp {src_file_bam} {dest_file_bam}\n")
        f.write(f"cp {src_file_bai} {dest_file_bai}\n")
