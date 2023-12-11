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
