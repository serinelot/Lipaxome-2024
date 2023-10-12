rule qc:
    """ Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_2.fastq.gz"
    
    output:
        qc_fq1_out = "data/qc/{id}_1_fastqc.html",
        qc_fq2_out = "data/qc/{id}_2_fastqc.html"
    
    params:
        out_dir = "data/qc"
    
    log:
        "logs/fastqc/{id}.log"
    
    threads:
        32
    
    conda:
        "../envs/fastqc.yml"
    
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}" 



rule trimming:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_1.fastq.gz"

    output:
        trimm_fq1 = "data/trimmed/{id}_1.fastq.gz",
        trimm_fq2 = "data/trimmed/{id}_2.fastq.gz",
        trimm_unpaired_fq1 = "data/trimmed/{id}_1.unpaired.fastq.gz",
        trimm_unpaired_fq2 = "data/trimmed/{id}_2.unpaired.fastq.gz"
    
    params:
        options = [
            "ILLUMINACLIP:data/Adapters-PE_NextSeq.fa:2:30:10",
            "LEADING:30", "TRAILING:30", "MINLEN:20"
        ]
    
    log:
        "logs/trimmomatic/{id}.log"
    
    threads:
        32
    
    conda:
        "../envs/trimmomatic.yml"
    
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq1} {input.fq2} "
        "{output.trimm_fq1} {output.trimm_unpaired_fq1}  "
        "{output.trimm_fq2} {output.trimm_unpaired_fq2} "
        "{params.options} "
        "&> {log}"



rule qc_trimm:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.trimming.output.trimm_fq1,
        trimm_fq2 = rules.trimming.output.trimm_fq2,
        trimm_unpaired_fq1 = rules.trimming.output.trimm_unpaired_fq1,
        trimm_unpaired_fq2 = rules.trimming.output.trimm_unpaired_fq2
    
    output:
        qc_trimm_fq1_out = "data/qc_trimmed/{id}_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_trimmed/{id}_2_fastqc.html"
    
    params:
        out_dir = "data/qc_trimmed"
    
    log:
        "logs/fastqc_trimmed/{id}.log"
    
    threads:
        32
    
    conda:
        "../envs/fastqc.yml"
    
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.trimm_fq1} {input.trimm_fq2} "
        "{input.trimm_unpaired_fq1} {input.trimm_unpaired_fq2} "
        "&> {log}" 