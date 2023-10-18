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
        fq2 = "data/fastq/{id}_2.fastq.gz"

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



rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = 'data/references/genome_fa/homo_sapiens_genome.fa',
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    output:
        chrNameLength = config['path']['chrNameLength']
    params:
        dir = config['path']['star_index']
    log:
        "logs/STAR/index.log"
    conda:
        "../envs/star.yml"
    threads:
        32
    shell:
        "mkdir -p {params.dir} && "
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99"
        "&> {log}"



rule star_alignReads:
    """ Generates a bam file using STAR """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trimming.output.trimm_fq1,
        fq2 = rules.trimming.output.trimm_fq2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/STAR/{id}/Log.final.out"
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/STAR/{id}.log"
    threads:
        32
    conda:
        "../envs/star.yml"
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.index} "
        "--readFilesIn {input.fq1} {input.fq2}  "
        "--runThreadN {threads} "
        "--readFilesCommand zcat "
        "--outReadsUnmapped Fastx "
        "--outFilterType BySJout "
        "--outStd Log "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.output_dir} "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 "
        "--winAnchorMultimapNmax 100 "
        "--limitBAMsortRAM 600000000000"
        "--alignEndsProtrude 5 ConcordantPair "
        "&> {log}"
