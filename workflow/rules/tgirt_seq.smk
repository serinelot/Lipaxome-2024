rule coco_ca:
    """ Generate corrected annotation from the gtf."""
    input:
        gtf = config["download"]["human_gtf"],
        coco_dir = rules.download_coco_git.output.git_coco_folder
    output:
        gtf_corrected = "data/references/gtf/hg38_Ensembl_V101_Scottlab_2020_correct_annotation.gtf"
    conda:
        "../envs/coco.yml"
    shell:
        "python3 git_repos/coco/bin/coco.py ca {input.gtf} -o {output.gtf_corrected}"



rule coco_cc:
    """Quantify the number of counts, counts per million (CPM) and transcript
        per million (TPM) for each gene using CoCo correct_count (cc)."""
    input:
        gtf = rules.coco_ca.output,
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        counts = "results/tgirt/coco_cc/{id}.tsv"
    threads:
        32
    params:
        coco_path = "git_repos/coco/bin"
    log:
        "logs/coco/coco_cc_{id}.log"
    conda:
        "../envs/coco.yml"
    shell:
        "python {params.coco_path}/coco.py cc "
        "--countType both "
        "--thread {threads} "
        "--strand 1 "
        "--paired "
        "{input.gtf} "
        "{input.bam} "
        "{output.counts} "
        "&> {log}"



rule merge_coco_cc_output:
    """ Merge CoCo correct count outputs into one count, cpm or tpm file (all
        samples merged inside one dataframe). This rule takes in input the
        output result directory of CoCo within the TGIRT-Seq pipeline."""
    input:
        counts = expand(rules.coco_cc.output.counts, id = id_list)
    output:
        merged_counts = "results/tgirt/coco_cc/merged_counts.tsv",
        merged_cpm = "results/tgirt/coco_cc/merged_cpm.tsv",
        merged_tpm = "results/tgirt/coco_cc/merged_tpm.tsv"
    conda:
        "../envs/coco.yml"
    script:
        "../scripts/merge_coco_cc_output.py"


rule deseq2_coco_mRNA:
    input:
        quant = expand("results/tgirt/coco_cc/{id}.tsv", id=id_list),
        samples = "data/design.tsv",
        comparisons = "data/comparisons.tsv",
        gene_id = "data/references/tx2gene.tsv"
    output:
        results = directory("results/dge/deseq2_coco_mRNA"),
        out_files = expand('results/dge/deseq2_coco_mRNA/{comp}.csv', comp = comparisons)
    params:
        quant_dir = "results/tgirt/coco_cc",
        filter_count_threshold = 10
    log:
        "logs/deseq2_coco_mRNA.log"
    message:
        "Perform differential expression analysis for various conditions."
    script:
        "../scripts/DESeq2_coco_tximport_mrna.R"



# rule coco_cb:
#     """ Create a bedgraph from the bam files """
#     input:
#         bam = rules.star_alignReads.output.bam,
#         chrNameLength = rules.star_index.output.chrNameLength
#     output:
#         unsorted_bedgraph = "results/tgirt/coco_cb/{id}_unsorted.bedgraph"
#     params:
#         pb = rules.install_pairedBamToBed12.params.pairedBamToBed12_bin,
#         coco_path = "git_repos/coco/bin"
#     conda:
#         "../envs/coco.yml"
#     threads:
#         28
#     shell:
#         "export PATH=$PWD/{params.pb}:$PATH && "
#         "python {params.coco_path}/coco.py cb "
#         "-u " # UCSC compatible (adds a chr and track info)
#         "-t {threads} "
#         "-c 2500000 " # Chunk size, default value
#         "{input.bam} "
#         "{output.unsorted_bedgraph} "
#         "{input.chrNameLength}"



# rule sort_bg:
#     """ Sort the new bedgraphs and change the chrM to chrMT to work with bedGraphToBigWig """
#     input:
#         unsorted_bedgraph = rules.coco_cb.output.unsorted_bedgraph
#     output:
#         sorted_bedgraph = "results/tgirt/coco_cb/{id}.bedgraph"
#     shell:
#         "sort -k1,1 -k2,2n {input.unsorted_bedgraph} "
#         "| sed 's/chrM/chrMT/g' > {output.sorted_bedgraph}"
