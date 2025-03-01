# Load config file
configfile: "conf.yaml"

# Define the sample list
SAMPLES = glob_wildcards("data/{sample}_R1.fastq.gz").sample
print(SAMPLES)
rule all:
    input:
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        expand("results/featurecounts/{sample}.counts.txt", sample=SAMPLES),
        expand("results/featurecounts/{sample}.counts.txt.summary", sample=SAMPLES)

# Rule: Quality Control with FastQC
rule fastqc:
    input:
        "data/{sample}_R1.fastq.gz",
        "data/{sample}_R2.fastq.gz"
    output:
        "qc/{sample}_R1_fastqc.html",
        "qc/{sample}_R2_fastqc.html"
    conda:
        "envs/rnaseq.yaml"
    shell:
        "fastqc {input} -o qc/"

# Rule: Adapter Trimming with Trim Galore
rule trim_galore:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        r1_trimmed="trimmed/{sample}_val_1.fq.gz",
        r2_trimmed="trimmed/{sample}_val_2.fq.gz"
    conda:
        "envs/rnaseq.yaml"
    shell:
        "trim_galore --paired {input.r1} {input.r2} --basename {wildcards.sample} --gzip -o trimmed/"

# Rule: Isoform based Quantification with Salmon
rule salmon_quant:
    input:
        r1="trimmed/{sample}_val_1.fq.gz",
        r2="trimmed/{sample}_val_2.fq.gz",
        index="references/salmon_index"  # Salmon index directory
    output:
        quant="results/salmon/{sample}/quant.sf"
    conda:
        "envs/rnaseq.yaml"
    shell:
        "{config[salmon_path]} quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} "
        "-o results/salmon/{wildcards.sample} --validateMappings"

# Rule: Alignment using STAR and Gene-level Quantification with featureCounts

rule align_reads:
    input:
        r1_trimmed="trimmed/{sample}_val_1.fq.gz",
        r2_trimmed="trimmed/{sample}_val_2.fq.gz"
    output:
        bam="aligned/{sample}.bam"
    conda:
        "envs/rnaseq.yaml"
    params:
        index="references/star_index"
    shell:
        """
        STAR --genomeDir {params.index} \
            --readFilesIn {input.r1_trimmed} {input.r2_trimmed} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix aligned/{wildcards.sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattrIHstart 1 \
            --outSAMprimaryFlag AllBestScore
        mv aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

rule featurecounts:
    input:
        bam="aligned/{sample}.bam",
        gtf="references/Homo_sapiens.GRCh38.110.gtf"
    output:
        counts="results/featurecounts/{sample}.counts.txt",
        summary="results/featurecounts/{sample}.counts.txt.summary"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        featureCounts -T 4 -p -t gene -g gene_id \
            -a {input.gtf} -o {output.counts} {input.bam}
        """
