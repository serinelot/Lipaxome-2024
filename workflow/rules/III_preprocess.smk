# workflow/rules/III_preprocess.smk

############################################################
# 1) FastQC avant trimming
############################################################
rule qc_pre_trim:
    """
    FastQC avant trimming : qualité des FASTQ bruts.
    """
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        html1 = "results/qc/pre_trim/{id}_1_fastqc.html",
        zip1  = "results/qc/pre_trim/{id}_1_fastqc.zip",
        html2 = "results/qc/pre_trim/{id}_2_fastqc.html",
        zip2  = "results/qc/pre_trim/{id}_2_fastqc.zip"
    params:
        outdir = "results/qc/pre_trim"
    log:
        "logs/fastqc/pre_trim/{id}.log"
    threads: 8
    message: "Running FastQC (pre-trim) for sample {wildcards.id}"
    conda:
        "../envs/fastqc.yml"
    shell:
        """
        fastqc \
          --outdir {params.outdir} \
          --format fastq \
          --threads {threads} \
          {input.fq1} {input.fq2} \
        &> {log}
        """

############################################################
# 2) Trimming (fastp)
############################################################
rule fastp_trim:
    """
    FastP : trimming des reads,
    outputs paired & unpaired FASTQ + rapports HTML/JSON.
    """
    input:
        fq1 = "data/fastq/{id}_1.fastq.gz",
        fq2 = "data/fastq/{id}_2.fastq.gz"
    output:
        paired1     = "results/preprocess/fastp/{id}_1_trimmed.fastq.gz",
        paired2     = "results/preprocess/fastp/{id}_2_trimmed.fastq.gz",
        unpaired1   = "results/preprocess/fastp/{id}_1_unpaired.fastq.gz",
        unpaired2   = "results/preprocess/fastp/{id}_2_unpaired.fastq.gz",
        html_report = "results/preprocess/fastp/{id}_fastp.html",
        json_report = "results/preprocess/fastp/{id}_fastp.json"
    params:
        outdir  = "results/preprocess/fastp",
        options = "--qualified_quality_phred 30 " \
                  "--length_required 20 " \
                  "--cut_window_size 1 " \
                  "--cut_mean_quality 30 " \
                  "--cut_front " \
                  "--cut_tail"
    log:
        "logs/fastp/{id}.log"
    threads: 8
    message: "Trimming reads for sample {wildcards.id}"
    conda:
        "../envs/fastp.yml"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastp \
          -i {input.fq1} -I {input.fq2} \
          -o {output.paired1} -O {output.paired2} \
          --unpaired1 {output.unpaired1} \
          --unpaired2 {output.unpaired2} \
          --thread {threads} \
          -h {output.html_report} \
          -j {output.json_report} \
          {params.options} \
        &> {log}
        """

############################################################
# 3) FastQC après trimming
############################################################
rule qc_post_trim:
    """
    FastQC après trimming : qualité des reads post‑fastp.
    """
    input:
        paired1   = rules.fastp_trim.output.paired1,
        paired2   = rules.fastp_trim.output.paired2,
        unpaired1 = rules.fastp_trim.output.unpaired1,
        unpaired2 = rules.fastp_trim.output.unpaired2
    output:
        html1 = "results/qc/post_trim/{id}_1_trimmed_fastqc.html",
        zip1  = "results/qc/post_trim/{id}_1_trimmed_fastqc.zip",
        html2 = "results/qc/post_trim/{id}_2_trimmed_fastqc.html",
        zip2  = "results/qc/post_trim/{id}_2_trimmed_fastqc.zip",
        html3 = "results/qc/post_trim/{id}_1_unpaired_fastqc.html",
        zip3  = "results/qc/post_trim/{id}_1_unpaired_fastqc.zip",
        html4 = "results/qc/post_trim/{id}_2_unpaired_fastqc.html",
        zip4  = "results/qc/post_trim/{id}_2_unpaired_fastqc.zip"
    params:
        outdir = "results/qc/post_trim"
    log:
        "logs/fastqc/post_trim/{id}.log"
    threads: 8
    message: "Running FastQC (post-trim) for sample {wildcards.id}"
    conda:
        "../envs/fastqc.yml"
    shell:
        r"""
        fastqc \
          --outdir {params.outdir} \
          --format fastq \
          --threads {threads} \
          {input.paired1} {input.paired2} {input.unpaired1} {input.unpaired2} \
        &> {log}
        """
