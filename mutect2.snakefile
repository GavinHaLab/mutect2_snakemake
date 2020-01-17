#mutect2.snakefile
#Anna Hoge
#Template made October 18th, 2019
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.2.4-foss-2016b-Python-3.6.6
ml java/jdk1.8.0_31
ml picard/2.18.29-Java

#command to run snakemake (remove -np at end when done validating):
snakemake -s mutect2.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"


rule all:
    input:
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.tbi", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}_f1r2.tar.gz", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.stats", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/mutect_merged.stats", tumors = config["samples"]),
        expand("results/{tumors}/read_orientation_model.tar.gz", tumors = config["samples"]),
        expand("results/{tumors}/pileup_summaries.table", tumors = config["samples"]),
        expand("results/{tumors}/segments.table", tumors = config["samples"]),
        expand("results/{tumors}/contamination.table", tumors = config["samples"]),
        expand("results/{tumors}/gathered_unfiltered.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/gathered_unfiltered.vcf.gz.tbi", tumors = config["samples"]),
        expand("results/{tumors}/filtered_all.vcf.gz", tumors = config["samples"]),
        expand("results/{tumors}/filtering_stats.tsv", tumors = config["samples"])


rule mutect2:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumors][0],
        normal_filepath = lambda wildcards: config["samples"][wildcards.tumors][2]
    output:
        vcf = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz"),
        tbi = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.tbi"),
        tar = protected("results/{tumors}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.stats")
    params:
        reference_genome = config["reference_genome"],
        mutect2_germline_resource = config["mutect2_germline_resource"],
        gatk = config["gatk"],
        panel_of_normals = config["panel_of_normals"],
        normals = lambda wildcards: config["samples"][wildcards.tumors][1]
    log:
        "logs/mutect2/{tumors}_{chromosomes}_mutect2.txt"
    shell:
        "({params.gatk} Mutect2 \
        -reference {params.reference_genome} \
        -input {input.tumor_filepath} \
        -input {input.normal_filepath} \
        -normal {params.normals} \
        -intervals {wildcards.chromosomes} \
        --germline-resource {params.mutect2_germline_resource} \
        --f1r2-tar-gz {output.tar} \
        --panel-of-normals {params.panel_of_normals} \
        -output {output.vcf}) 2> {log}"


rule merge_mutect_stats:
    input:
        chr1_stats = "results/{tumors}/unfiltered_chr1.vcf.gz.stats",
        chr2_stats = "results/{tumors}/unfiltered_chr2.vcf.gz.stats",
        chr3_stats = "results/{tumors}/unfiltered_chr3.vcf.gz.stats",
        chr4_stats = "results/{tumors}/unfiltered_chr4.vcf.gz.stats",
        chr5_stats = "results/{tumors}/unfiltered_chr5.vcf.gz.stats",
        chr6_stats = "results/{tumors}/unfiltered_chr6.vcf.gz.stats",
        chr7_stats = "results/{tumors}/unfiltered_chr7.vcf.gz.stats",
        chr8_stats = "results/{tumors}/unfiltered_chr8.vcf.gz.stats",
        chr9_stats = "results/{tumors}/unfiltered_chr9.vcf.gz.stats",
        chr10_stats = "results/{tumors}/unfiltered_chr10.vcf.gz.stats",
        chr11_stats = "results/{tumors}/unfiltered_chr11.vcf.gz.stats",
        chr12_stats = "results/{tumors}/unfiltered_chr12.vcf.gz.stats",
        chr13_stats = "results/{tumors}/unfiltered_chr13.vcf.gz.stats",
        chr14_stats = "results/{tumors}/unfiltered_chr14.vcf.gz.stats",
        chr15_stats = "results/{tumors}/unfiltered_chr15.vcf.gz.stats",
        chr16_stats = "results/{tumors}/unfiltered_chr16.vcf.gz.stats",
        chr17_stats = "results/{tumors}/unfiltered_chr17.vcf.gz.stats",
        chr18_stats = "results/{tumors}/unfiltered_chr18.vcf.gz.stats",
        chr19_stats = "results/{tumors}/unfiltered_chr19.vcf.gz.stats",
        chr20_stats = "results/{tumors}/unfiltered_chr20.vcf.gz.stats",
        chr21_stats = "results/{tumors}/unfiltered_chr21.vcf.gz.stats",
        chr22_stats = "results/{tumors}/unfiltered_chr22.vcf.gz.stats",
        chrX_stats = "results/{tumors}/unfiltered_chrX.vcf.gz.stats",
        chrY_stats = "results/{tumors}/unfiltered_chrY.vcf.gz.stats"
    output:
        protected("results/{tumors}/mutect_merged.stats")
    params:
        gatk = config["gatk"]
    log:
        "logs/merge_mutect_stats/{tumors}_merge_mutect_stats.txt"
    shell:
        "({params.gatk} MergeMutectStats \
        -stats {input.chr1_stats} \
        -stats {input.chr2_stats} \
        -stats {input.chr3_stats} \
        -stats {input.chr4_stats} \
        -stats {input.chr5_stats} \
        -stats {input.chr6_stats} \
        -stats {input.chr7_stats} \
        -stats {input.chr8_stats} \
        -stats {input.chr9_stats} \
        -stats {input.chr10_stats} \
        -stats {input.chr11_stats} \
        -stats {input.chr12_stats} \
        -stats {input.chr13_stats} \
        -stats {input.chr14_stats} \
        -stats {input.chr15_stats} \
        -stats {input.chr16_stats} \
        -stats {input.chr17_stats} \
        -stats {input.chr18_stats} \
        -stats {input.chr19_stats} \
        -stats {input.chr20_stats} \
        -stats {input.chr21_stats} \
        -stats {input.chr22_stats} \
        -stats {input.chrX_stats} \
        -stats {input.chrY_stats} \
        -O {output}) 2> {log}"

rule learn_read_orientation_model:
    input:
        chr1_tar = "results/{tumors}/unfiltered_chr1_f1r2.tar.gz",
        chr2_tar = "results/{tumors}/unfiltered_chr2_f1r2.tar.gz",
        chr3_tar = "results/{tumors}/unfiltered_chr3_f1r2.tar.gz",
        chr4_tar = "results/{tumors}/unfiltered_chr4_f1r2.tar.gz",
        chr5_tar = "results/{tumors}/unfiltered_chr5_f1r2.tar.gz",
        chr6_tar = "results/{tumors}/unfiltered_chr6_f1r2.tar.gz",
        chr7_tar = "results/{tumors}/unfiltered_chr7_f1r2.tar.gz",
        chr8_tar = "results/{tumors}/unfiltered_chr8_f1r2.tar.gz",
        chr9_tar = "results/{tumors}/unfiltered_chr9_f1r2.tar.gz",
        chr10_tar = "results/{tumors}/unfiltered_chr10_f1r2.tar.gz",
        chr11_tar = "results/{tumors}/unfiltered_chr11_f1r2.tar.gz",
        chr12_tar = "results/{tumors}/unfiltered_chr12_f1r2.tar.gz",
        chr13_tar = "results/{tumors}/unfiltered_chr13_f1r2.tar.gz",
        chr14_tar = "results/{tumors}/unfiltered_chr14_f1r2.tar.gz",
        chr15_tar = "results/{tumors}/unfiltered_chr15_f1r2.tar.gz",
        chr16_tar = "results/{tumors}/unfiltered_chr16_f1r2.tar.gz",
        chr17_tar = "results/{tumors}/unfiltered_chr17_f1r2.tar.gz",
        chr18_tar = "results/{tumors}/unfiltered_chr18_f1r2.tar.gz",
        chr19_tar = "results/{tumors}/unfiltered_chr19_f1r2.tar.gz",
        chr20_tar = "results/{tumors}/unfiltered_chr20_f1r2.tar.gz",
        chr21_tar = "results/{tumors}/unfiltered_chr21_f1r2.tar.gz",
        chr22_tar = "results/{tumors}/unfiltered_chr22_f1r2.tar.gz",
        chrX_tar = "results/{tumors}/unfiltered_chrX_f1r2.tar.gz",
        chrY_tar = "results/{tumors}/unfiltered_chrY_f1r2.tar.gz"
    output:
        protected("results/{tumors}/read_orientation_model.tar.gz")
    params:
        gatk = config["gatk"]
    log:
        "logs/learn_read_orientation_model/{tumors}_learn_read_orientation_model.txt"
    shell:
        "({params.gatk} LearnReadOrientationModel \
        -I {input.chr1_tar} \
        -I {input.chr2_tar} \
        -I {input.chr3_tar} \
        -I {input.chr4_tar} \
        -I {input.chr5_tar} \
        -I {input.chr6_tar} \
        -I {input.chr7_tar} \
        -I {input.chr8_tar} \
        -I {input.chr9_tar} \
        -I {input.chr10_tar} \
        -I {input.chr11_tar} \
        -I {input.chr12_tar} \
        -I {input.chr13_tar} \
        -I {input.chr14_tar} \
        -I {input.chr15_tar} \
        -I {input.chr16_tar} \
        -I {input.chr17_tar} \
        -I {input.chr18_tar} \
        -I {input.chr19_tar} \
        -I {input.chr20_tar} \
        -I {input.chr21_tar} \
        -I {input.chr22_tar} \
        -I {input.chrX_tar} \
        -I {input.chrY_tar} \
        -O {output}) 2> {log}"


rule get_pileup_summaries:
    input:
        lambda wildcards: config["samples"][wildcards.tumors][0]
    output:
        protected("results/{tumors}/pileup_summaries.table")
    params:
        gatk = config["gatk"],
        known_polymorphic_sites = config["known_polymorphic_sites"]
    log:
        "logs/get_pileup_summaries/{tumors}_get_pileup_summaries.txt"
    shell:
        "({params.gatk} GetPileupSummaries \
        -I {input} \
        -V {params.known_polymorphic_sites} \
        -L {params.known_polymorphic_sites} \
        -O {output}) 2> {log}"


rule calculate_contamination:
    input:
        "results/{tumors}/pileup_summaries.table"
    output:
        segments_table = protected("results/{tumors}/segments.table"),
        contamination_table = protected("results/{tumors}/contamination.table")
    params:
        gatk = config["gatk"]
    log:
        "logs/calculate_contamination/{tumors}_calculate_contamination.txt"
    shell:
        "({params.gatk} CalculateContamination \
        -I {input} \
        -tumor-segmentation {output.segments_table} \
        -O {output.contamination_table}) 2> {log}"


rule gather_mutect_calls:
    input:
        chr1_calls = "results/{tumors}/unfiltered_chr1.vcf.gz",
        chr2_calls = "results/{tumors}/unfiltered_chr2.vcf.gz",
        chr3_calls = "results/{tumors}/unfiltered_chr3.vcf.gz",
        chr4_calls = "results/{tumors}/unfiltered_chr4.vcf.gz",
        chr5_calls = "results/{tumors}/unfiltered_chr5.vcf.gz",
        chr6_calls = "results/{tumors}/unfiltered_chr6.vcf.gz",
        chr7_calls = "results/{tumors}/unfiltered_chr7.vcf.gz",
        chr8_calls = "results/{tumors}/unfiltered_chr8.vcf.gz",
        chr9_calls = "results/{tumors}/unfiltered_chr9.vcf.gz",
        chr10_calls = "results/{tumors}/unfiltered_chr10.vcf.gz",
        chr11_calls = "results/{tumors}/unfiltered_chr11.vcf.gz",
        chr12_calls = "results/{tumors}/unfiltered_chr12.vcf.gz",
        chr13_calls = "results/{tumors}/unfiltered_chr13.vcf.gz",
        chr14_calls = "results/{tumors}/unfiltered_chr14.vcf.gz",
        chr15_calls = "results/{tumors}/unfiltered_chr15.vcf.gz",
        chr16_calls = "results/{tumors}/unfiltered_chr16.vcf.gz",
        chr17_calls = "results/{tumors}/unfiltered_chr17.vcf.gz",
        chr18_calls = "results/{tumors}/unfiltered_chr18.vcf.gz",
        chr19_calls = "results/{tumors}/unfiltered_chr19.vcf.gz",
        chr20_calls = "results/{tumors}/unfiltered_chr20.vcf.gz",
        chr21_calls = "results/{tumors}/unfiltered_chr21.vcf.gz",
        chr22_calls = "results/{tumors}/unfiltered_chr22.vcf.gz",
        chrX_calls = "results/{tumors}/unfiltered_chrX.vcf.gz",
        chrY_calls = "results/{tumors}/unfiltered_chrY.vcf.gz"
    output:
        protected("results/{tumors}/gathered_unfiltered.vcf.gz")
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/gather_mutect_calls/{tumors}_gather_mutect_calls.txt"
    shell:
        "({params.java} -jar {params.picard_jar} GatherVcfs \
        I={input.chr1_calls} \
        I={input.chr2_calls} \
        I={input.chr3_calls} \
        I={input.chr4_calls} \
        I={input.chr5_calls} \
        I={input.chr6_calls} \
        I={input.chr7_calls} \
        I={input.chr8_calls} \
        I={input.chr9_calls} \
        I={input.chr10_calls} \
        I={input.chr11_calls} \
        I={input.chr12_calls} \
        I={input.chr13_calls} \
        I={input.chr14_calls} \
        I={input.chr15_calls} \
        I={input.chr16_calls} \
        I={input.chr17_calls} \
        I={input.chr18_calls} \
        I={input.chr19_calls} \
        I={input.chr20_calls} \
        I={input.chr21_calls} \
        I={input.chr22_calls} \
        I={input.chrX_calls} \
        I={input.chrY_calls} \
        O={output}) 2> {log}"


rule create_vcf_index:
    input:
        "results/{tumors}/gathered_unfiltered.vcf.gz"
    output:
        protected("results/{tumors}/gathered_unfiltered.vcf.gz.tbi")
    log:
        "logs/create_vcf_index/{tumors}_create_vcf_index.txt"
    shell:
        "(tabix -p vcf {input}) 2> {log}"


rule filter_mutect_calls:
    input:
        unfiltered_vcf = "results/{tumors}/gathered_unfiltered.vcf.gz",
        vcf_index = "results/{tumors}/gathered_unfiltered.vcf.gz.tbi",
        segments_table = "results/{tumors}/segments.table",
        contamination_table = "results/{tumors}/contamination.table",
        read_orientation_model = "results/{tumors}/read_orientation_model.tar.gz",
        mutect_stats = "results/{tumors}/mutect_merged.stats"
    output:
        filtered_vcf = protected("results/{tumors}/filtered_all.vcf.gz"),
        filtering_stats = protected("results/{tumors}/filtering_stats.tsv")
    params:
        gatk = config["gatk"],
        reference_genome = config["reference_genome"]
    log:
        "logs/filter_mutect_calls/{tumors}_filter_mutect_calls.txt"
    shell:
        "({params.gatk} FilterMutectCalls \
        -R {params.reference_genome} \
        -V {input.unfiltered_vcf} \
        --tumor-segmentation {input.segments_table} \
        --contamination-table {input.contamination_table} \
        --ob-priors {input.read_orientation_model} \
        --stats {input.mutect_stats} \
        --filtering-stats {output.filtering_stats} \
        -O {output.filtered_vcf}) 2> {log}"
