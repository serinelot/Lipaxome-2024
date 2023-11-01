import os

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



rule download_human_sample_fastq:
    """Download fastq of human samples from GEO. This step
        might take long, depending on your downloading speed. Other samples were
        manually added on the cluster. To automate the download of ALL samples,
        simply add their SRR id in the config file (use SRA explorer to find the URL)"""
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



rule install_pairedBamToBed12:
    output:
        directory("scripts/pairedBamToBed12/")
    params:
        pairedBamToBed12_bin = 'scripts/pairedBamToBed12/bin'
    conda:
        "../envs/coco.yml"
    shell:
        'mkdir -p scripts && cd scripts && pwd && '
        'git clone https://github.com/Population-Transcriptomics/pairedBamToBed12 && '
        'cd pairedBamToBed12 && '
        'make '




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
