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
