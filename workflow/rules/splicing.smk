rule primary_alignments:
    # More robust results with multimapped reads
    input:
        rules.star_alignReads.output.bam
    output:
        "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam"
    log:
        "logs/STAR/{id}_primary.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Keep primary alignments only for {wildcards.id}."
    shell:
        "samtools view -b -F 256 -o {output} {input} "
        "&> {log}"



rule bam_index:
    input:
        rules.primary_alignments.output
    output:
        "results/splicing/star/{id}/{id}_Aligned.sortedByCoord.out.primary.bam.bai"
    log:
        "logs/star/{id}_primary_index.log"
    conda:
        "../envs/genomecov.yml"
    message:
        "Create a BAI index for {wildcards.id}."
    shell:
        "samtools index {input} "
        "&> {log}"



rule genomecov:
    input:
        rules.primary_alignments.output
    output:
        "results/splicing/genomecov/{id}.bedgraph"
    conda:
        "../envs/genomecov.yml"
    message:
        "Report {wildcards.id} genome coverage in BEDGRAPH format."
    shell:
        "bedtools genomecov -bg -split -ibam {input} | sort -k1,1 -k2,2n > {output}" 
