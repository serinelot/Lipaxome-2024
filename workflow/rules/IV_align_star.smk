rule star_index:
    """
    Génère l’index genome pour STAR à partir du FASTA et du GTF.
    """
    input:
        fasta = rules.download_human_genome.output.genome,
        gtf   = config["download"]["human_gtf"]
    output:
        chrNameLength = config["path"]["chrNameLength"]
    params:
        index_dir = config["path"]["star_index"]
    log:
        "logs/STAR/index.log"
    threads: 32
    message: "Building STAR index in {params.index_dir}"
    conda:
        "../envs/star.yml"
    shell:
        r"""
        mkdir -p {params.index_dir}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {params.index_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 99 \
        &> {log}
        """
        
rule star_alignReads:
    """
    Génère un BAM trié et un Log.final.out via STAR.
    Utilisé pour fournir les alignements aux étapes de quantification.
    """
    input:
        idx = rules.star_index.output.chrNameLength,
        fq1 = rules.fastp_trim.output.paired1,
        fq2 = rules.fastp_trim.output.paired2
    output:
        bam      = "results/alignment/star/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/alignment/star/{id}/Log.final.out"
    params:
        index      = config["path"]["star_index"],
        output_dir = "results/alignment/star/{id}/"
    log:
        "logs/STAR/align_{id}.log"
    threads: 32
    message: "Aligning sample {wildcards.id} with STAR"
    conda:
        "../envs/star.yml"
    shell:
        r"""
        mkdir -p {params.output_dir}
        STAR --runMode alignReads \
             --genomeDir {params.index} \
             --readFilesIn {input.fq1} {input.fq2} \
             --runThreadN {threads} \
             --readFilesCommand zcat \
             --outReadsUnmapped Fastx \
             --outFilterType BySJout \
             --outStd Log \
             --outSAMunmapped None \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.output_dir} \
             --outFilterScoreMinOverLread 0.3 \
             --outFilterMatchNminOverLread 0.3 \
             --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 100 \
             --limitBAMsortRAM 600000000000 \
             --alignEndsProtrude 5 ConcordantPair \
        &> {log}
        """
