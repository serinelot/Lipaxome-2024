rule build_transcriptome:
    input:
        genome = rules.download_human_genome.output.genome,
        gtf = config['download']['human_gtf']   
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
        fq1 = rules.fastp.output.fastq1,
        fq2 = rules.fastp.output.fastq2
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
        gtf = config['download']['human_gtf']
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
        gtf = config['download']['human_gtf']
    output:
        pc_gtf = "data/references/gtf/Homo_sapiens.GRCh38.110_snoRNAs_tRNAs.protein_coding.gtf"
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


rule deseq2_mRNA:
    input:
        quant = expand("results/dge/kallisto/{id}/abundance.tsv", id=id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = "data/references/tx2gene.tsv"
    output:
        results = directory("results/dge/deseq2_kallisto_mRNA"),
        out_files = expand('results/dge/deseq2_kallisto_mRNA/{comp}.csv', comp = comparisons)
    params:
        kallisto_dir = "results/dge/kallisto",
        filter_count_threshold = 10
    log:
        "logs/deseq2_mRNA.log"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_kallisto_tximport_mrna.R"


rule deseq2_mrna_stat:
    input:
        deseq2 = rules.deseq2_mRNA.output.out_files
    output:
        deseq2_stat = expand('results/dge/deseq2_kallisto_mRNA/{comp}_deseq2_stat.csv', comp = comparisons), 
        volcano_plot_total = expand('results/dge/deseq2_kallisto_mRNA/{comp}_volcano_plot_total.png', comp = comparisons),
        volcano_plot_zoom = expand('results/dge/deseq2_kallisto_mRNA/{comp}_volcano_plot_zoom.png', comp = comparisons),
        go_enrichment_total = expand('results/dge/deseq2_kallisto_mRNA/{comp}_go_enrichment_total.png', comp = comparisons),
        go_enrichment_up = expand('results/dge/deseq2_kallisto_mRNA/{comp}_go_enrichment_up.png', comp = comparisons),
        go_enrichment_down = expand('results/dge/deseq2_kallisto_mRNA/{comp}_go_enrichment_down.png', comp = comparisons)
    log:
        "logs/deseq2_mRNA_stat.log"
    script:
        "../scripts/deseq2_stat_kallisto_mRNA.R"


rule deseq2_mRNA_transcript:
    input:
        quant = expand("results/dge/kallisto/{id}/abundance.tsv", id=id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = "data/references/tx2gene.tsv"
    output:
        results = directory("results/dge/deseq2_kallisto_mRNA_transcript"),
        out_files = expand('results/dge/deseq2_kallisto_mRNA_transcript/{comp}.csv', comp = comparisons)
    params:
        kallisto_dir = "results/dge/kallisto",
        filter_count_threshold = 10
    log:
        "logs/deseq2_mRNA_trancript.log"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_kallisto_tximport_mrna_transcript.R"

# rule volcano_plot:
#     input:
#         DE_outdir = rules.deseq2.output.results,
#         DE_output = "results/dge/deseq2_kallisto_mRNA/{comp}.csv",
#         filtered_genes = rules.merge_kallisto_quant.output.tpm,
#         gtf = rules.filter_gtf_pc_genes.output.pc_gtf
#     output:
#         volcano = "results/dge/volcano_plot/deseq2_kallisto_mrna_{comp}.svg",
#         up_genes = "results/dge/volcano_plot/deseq2_kallisto_mrna_{comp}_up_genes.tsv",
#         down_genes = "results/dge/volcano_plot/deseq2_kallisto_mrna_{comp}_down_genes.tsv"
#     params:
#         pval_threshold = 0.05
#     log:
#         "logs/volcano_plot/{comp}.log"
#     conda:
#         "../envs/python.yml"
#     message:
#         "Create a volcano plot using deseq2 output for each comparison in comparisons.tsv."
#     script:
#         "../scripts/volcano_plot.py"


# rule GO_analysis_upregulated_genes:
#     input:
#         genes = rules.volcano_plot.output.up_genes,
#         go_obo = rules.go_basic_obo.output.go_obo,
#         go_gaf = rules.go_annotation_gaf.output.go_gaf
#     output:
#         bar_chart = "results/dge/GO/GO_{comp}_upregulated_genes.svg"
#     log:
#         "logs/GO/{comp}_upregulated_genes.log"
#     conda: 
#         "../envs/GO.yml"
#     message:
#         "GO analysis of upregulated genes in {wildcards.comp} represented as a bar chart."
#     script:
#         "../scripts/GO_bar_charts.py"
