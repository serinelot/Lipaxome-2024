rule download_human_genome:
    """
    Télécharger le génome de référence humain (fichier fasta)
    depuis les serveurs FTP d’ENSEMBL, le décompresser
    puis le déplacer à son emplacement final.
    """
    output:
        genome = "data/references/genome_fa/Homo_sapiens.GRCh38.115.dna.primary_assembly.fa"
    params:
        link = config["download"]["human_genome_fa"]
    log:
        "logs/download/genome_fa.log"
    message: "Downloading and preparing human reference genome"
    shell:
        r"""
        mkdir -p data/references/genome_fa
        wget -O data/references/genome_fa/temp.fa.gz {params.link} \
          &> {log}
        gunzip -f data/references/genome_fa/temp.fa.gz \
          &>> {log}
        mv data/references/genome_fa/temp.fa {output.genome}
        """

rule fai_to_chromsizes:
    """
    Génère le fichier .chrom.sizes à partir de l'index fasta (.fai).
    """
    input:
        fai = "data/references/genome_fa/Homo_sapiens.GRCh38.115.dna.primary_assembly.fa.fai"
    output:
        chrom_sizes = "data/references/genome_fa/GRCh38.115.chrom.sizes"
    message:
        "Extraction des tailles de chromosomes depuis le .fai"
    shell:
        """
        cut -f1,2 {input.fai} > {output.chrom_sizes}
        """


rule clone_coco:
    """
    Cloner le dépôt COCO pour la quantification des reads.
    """
    output:
        repo_dir = directory("git_repos/coco")
    params:
        url = config["path"]["coco_git_link"]
    log:
        "logs/download/coco_clone.log"
    message:
        "Cloning COCO repository"
    shell:
        r"""
        mkdir -p git_repos
        git clone {params.url} {output.repo_dir} &> {log}
        """

rule tx2gene:
    """
    Générer localement tx2gene.tsv à partir du GTF.
    Utile pour regrouper les abondances transcriptales en abondances géniques.
    """
    input:
        gtf = config["download"]["human_gtf"]
    output:
        tsv = "data/references/tx2gene.tsv"
    log:
        "logs/download/tx2gene.log"
    message:
        "Building transcript-to-gene map"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/tx2gene.py"


rule build_tx2gene_all:
    """
    Générer tx2gene_all.tsv (tous biotypes) depuis le GTF complet
    Utile pour sommer isoformes coding + non-coding.
    """
    input:
        gtf = config["download"]["human_gtf"]
    output:
        tx2gene_all = "data/references/tx2gene_all.tsv"
    log:
        "logs/references/tx2gene_all.log"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/tx2gene_all.py"
