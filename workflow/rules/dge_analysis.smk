rule build_transcriptome:
    input:
        genome = config["path"]["human_genome_fa"],
        gtf = config["path"]["human_gtf"]    
    output:
        config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yml"
    message:
        "Build a reference transcriptome using gffread."
    log:
        "logs/build_transcriptome/build_transcriptome.log"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output}"



rule kallisto_index:
    input:
        rules.build_transcriptome.output
    output:
        "data/references/kallisto.idx"
    params:
        31
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/index.log"
    message:
        "Builds an index from the FASTA file."
    shell:
        "kallisto index "
        "--index={output} "
        "--kmer-size={params} "
        "{input} "
        "&> {log}"



rule kallisto_quant:
    input:
        idx = rules.kallisto_index.output,
        fq1 = rules.trimming.output.trimm_fq1,
        fq2 = rules.trimming.output.trimm_fq2
    output:
        "results/dge/kallisto/{id}/abundance.tsv"
    params:
        bootstrap = "50",
        outdir = "results/dge/kallisto/{id}"
    threads:
        1
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/{id}.log"
    message:
        "Perform pseudoalignment and quantify transcript abundance for {wildcards.id}."
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"



rule tx2gene:
    input:
        gtf = config["path"]["human_gtf"]
    output:
        tsv = "data/references/tx2gene.tsv"
    conda:
        "../envs/python.yml"
    message:
        "Convert transcript IDs to gene IDs."
    script:
        "../scripts/tx2gene.py"



rule filter_gtf_pc_genes:
    input:
        gtf = config["path"]["human_gtf"]
    output:
        pc_gtf = "data/references/gtf/hg38_Ensembl_V101_Scottlab_2020.protein_coding.gtf"
    log:
        "logs/kallisto/filter_gtf_pc_genes.log"
    message:
        "Extract protein coding genes from the genome annotation file."
    shell:
        "grep \'protein_coding\' {input} > {output}"



rule merge_kallisto_quant:
    input:
        quant = expand(rules.kallisto_quant.output, id = id_list),
        tx2gene = rules.tx2gene.output.tsv,
        gtf = rules.filter_gtf_pc_genes.output.pc_gtf
    output:
        tpm = "results/dge/kallisto/tpm.tsv"
    conda:
        "../envs/python.yml"
    log:
        "logs/kallisto/merge_kallisto_quant.log"
    message:
        "Merge kallisto quantification results into one dataframe for further analysis."
    script:
        "../scripts/merge_kallisto_quant.py"




rule extract_cdn_values:
    output:
        cdn1_file = "data/cdn1.txt",
        cdn2_file = "data/cdn2.txt"
    run:
        with open("data/comparisons.tsv") as f, open(output.cdn1_file, 'w') as f1, open(output.cdn2_file, 'w') as f2:
            lines = f.read().strip().split('\n')
            cdn1 = lines[1].split()[0]  # Assuming the first line contains column names
            cdn2 = lines[1].split()[1]
            f1.write(cdn1)
            f2.write(cdn2)




rule deseq2:
    input:
        quant = expand("results/dge/kallisto/{id}/abundance.tsv", id = id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = "data/references/tx2gene.tsv"
    output:
        results = directory("results/dge/deseq2"),
        deseq2_output = "results/dge/deseq2/FXS-Control.csv"
    params:
        kallisto_dir = "results/dge/kallisto"
    log:
        "logs/deseq2.log"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_kallisto_tximport.R"
