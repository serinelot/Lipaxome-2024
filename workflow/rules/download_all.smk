rule download_human_sample_fastq:
    """Téléchargement des données de RNA-Seq du projet Lipaxome"""
    output:
        "data/references/geo_download.txt"
    params:
        sample_list = "data/references/sra_id_wget.txt"
    shell:
        "mkdir -p data/fastq/ && "
        "cd data/fastq/ && "
        "wget -i ../references/sra_id_wget.txt && "
        "cd ../references/ && "
        "touch geo_download.txt"


# rule download_human_gtf:
#     """ Download gtf of human genome from Ensembl """
#     output:
#         gtf = 'data/references/gtf/homo_sapiens.gtf'
#     params:
#         link = config['download']['human_gtf']
#     shell:
#         "wget -O temp.gz {params.link} && "
#         "gunzip temp.gz && "
#         "mv temp {output.gtf}"


rule download_human_gff3:
    """ Download gff3 of human genome from Ensembl """
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
        link = config['download']['human_genome_fa']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"


# rule install_pairedBamToBed12:
#     output:
#         directory("scripts/pairedBamToBed12/")
#     params:
#         pairedBamToBed12_bin = 'scripts/pairedBamToBed12/bin'
#     conda:
#         "../envs/coco.yml"
#     shell:
#         'mkdir -p scripts && cd scripts && pwd && '
#         'git clone https://github.com/Population-Transcriptomics/pairedBamToBed12 && '
#         'cd pairedBamToBed12 && '
#         'make '


rule download_coco_git:
    """Download git repository of CoCo."""
    output:
        git_coco_folder = directory('git_repos/coco')
    params:
        git_coco_link = config['path']['coco_git_link']
    conda:
        '../envs/git.yml'
    shell:
        'mkdir -p {output.git_coco_folder} '
        '&& git clone {params.git_coco_link} {output.git_coco_folder}'


rule go_basic_obo:
    """ Download go-basic.obo file from Gene Ontology """
    output:
        go_obo = 'data/references/GO/go-basic.obo'
    params:
        link = config['download']['go_basic_obo']
    shell:
        "wget -O temp {params.link} && "
        "mv temp {output.go_obo}"


rule go_annotation_gaf:
    """ Download goa_human_rna.gaf file from Gene Ontology """
    output:
        go_gaf = 'data/references/GO/goa_human_rna.gaf'
    params:
        link = config['download']['go_annotation_gaf']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.go_gaf}"
