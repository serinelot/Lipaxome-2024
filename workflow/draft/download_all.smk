###############################################################################################################################
#################### À UTILISER SI ON DOIT TÉLÉCHARGER CES FICHIERS (ET NON UTILISER LES FICHIERS DU LABO SCOTT) ##############
###############################################################################################################################

rule download_human_gtf:
    """ Download gtf of human genome from Zenodo"""
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"



rule download_human_gff3:
    """ Download gff3 of human genome """
    output:
        gff3 = 'data/references/gff3/homo_sapiens.gff3'
    params:
        link = config['download']['human_gff3']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gff3}"



rule download_human_genome:
    """Download the reference genome (fasta file) of human
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/homo_sapiens_genome.fa'
    params:
        link = "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"
