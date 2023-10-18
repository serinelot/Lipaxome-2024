rule feature_counts:
    input:
        bam = rules.star_alignReads.output.bam, 
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    output:
        counts = "results/featureCounts/{id}_gene_counts.tsv"
    log:
        "logs/featureCounts/{id}.log"
    threads:
        32
    conda:
        "../envs/featureCounts.yml"
    shell:
        "featureCounts "
        "-a {input.gtf} "
        "-g gene_id "
        "-M " 
        "-o {output.counts} "
        "-p "
        "-s 1 "
        "-T {threads} "
        "{input.bam} "
