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


rule trimmomatic:
    """ Trims the FASTQ files using Trimmomatic """
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        trimm_fq1 = "data/trimmomatic/{id}_1.fastq.gz",
        trimm_fq2 = "data/trimmomatic/{id}_2.fastq.gz",
        trimm_unpaired_fq1 = "data/trimmomatic/{id}_1.unpaired.fastq.gz",
        trimm_unpaired_fq2 = "data/trimmomatic/{id}_2.unpaired.fastq.gz"
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


rule qc_trimmomatic:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.trimmomatic.output.trimm_fq1,
        trimm_fq2 = rules.trimmomatic.output.trimm_fq2,
        trimm_unpaired_fq1 = rules.trimmomatic.output.trimm_unpaired_fq1,
        trimm_unpaired_fq2 = rules.trimmomatic.output.trimm_unpaired_fq2
    output:
        qc_trimm_fq1_out = "data/qc_trimmomatic/{id}_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_trimmomatic/{id}_2_fastqc.html"
    params:
        out_dir = "data/qc_trimmomatic"
    log:
        "logs/qc_trimmomatic/{id}.log"
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


rule fastp:
    """ Trim reads from fastq files using fastp."""
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        fastq1 = "data/fastp/{id}_1.fastq.gz",
        fastq2 = "data/fastp/{id}_2.fastq.gz",
        unpaired_fastq1 = "data/fastp/{id}_1.unpaired.fastq.gz",
        unpaired_fastq2 = "data/fastp/{id}_2.unpaired.fastq.gz",
        html_report = "data/fastp/{id}_fastp.html",
        json_report = "data/fastp/{id}_fastp.json"
    threads:
        8
    params:
        options = ["--qualified_quality_phred 30", "--length_required 20",
                "--cut_window_size 1", "--cut_mean_quality 30", "--cut_front",
                "--cut_tail"]
    log:
        "logs/fastp/{id}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        "fastp -i {input.fq1} -I {input.fq2} "
        "-o {output.fastq1} -O {output.fastq2} "
        "--unpaired1 {output.unpaired_fastq1} "
        "--unpaired2 {output.unpaired_fastq2} "
        "--thread {threads} "
        "-h {output.html_report} "
        "-j {output.json_report} "
        "{params.options} "
        "&> {log}"


rule qc_fastp:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.fastp.output.fastq1,
        trimm_fq2 = rules.fastp.output.fastq2,
        trimm_unpaired_fq1 = rules.fastp.output.unpaired_fastq1,
        trimm_unpaired_fq2 = rules.fastp.output.unpaired_fastq2
    output:
        qc_trimm_fq1_out = "data/qc_fastp/{id}_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_fastp/{id}_2_fastqc.html"
    params:
        out_dir = "data/qc_fastp"
    log:
        "logs/qc_fastp/{id}.log"
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
        fasta = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf
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
        fq1 = rules.fastp.output.fastq1,
        fq2 = rules.fastp.output.fastq2
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
