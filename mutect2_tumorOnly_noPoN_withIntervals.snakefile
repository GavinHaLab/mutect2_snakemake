#mutect2_tumorOnly_noPoN_withIntervals.snakefile
#Modified version: NO Panel of Normals, WITH Intervals file, WITH ANNOVAR
"""

ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Java/11.0.2
ml picard/2.21.6-Java-11
ml GATK/4.1.8.1-GCCcore-8.3.0-Java-11
ml tabix/0.2.6-GCCcore-8.3.0
ml annovar

#command to run snakemake (remove -np at end when done validating):
snakemake -s mutect2_tumorOnly_noPoN_withIntervals.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/samples_tumorOnly.yaml"
configfile: "config/config_tumorOnly_with_Intervals.yaml"

rule all:
    input:
        expand("results_tumorOnly/{tumors}/unfiltered.vcf.gz.stats", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/read_orientation_model.tar.gz", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/pileup_summaries.table", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/segments.table", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/contamination.table", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/filtered_all.vcf.gz", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/filtering_stats.tsv", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/pass_variants.vcf.gz", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/pass_snvs.vcf.gz", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/pass_variants.hg38_multianno.vcf", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/pass_snvs.hg38_multianno.vcf", tumors=config["samples"]),
        expand("results_tumorOnly/{tumors}/filtered_all.hg38_multianno.vcf", tumors=config["samples"])

rule mutect2_tumor_only:
    input:
        tumor=lambda wc: config["samples"][wc.tumors]
    output:
        vcf=temp("results_tumorOnly/{tumors}/unfiltered.vcf.gz"),
        tbi=temp("results_tumorOnly/{tumors}/unfiltered.vcf.gz.tbi"),
        tar=temp("results_tumorOnly/{tumors}/unfiltered.f1r2.tar.gz"),
        stats=protected("results_tumorOnly/{tumors}/unfiltered.vcf.gz.stats")
    params:
        reference=config["reference_genome"],
        germline=config["mutect2_germline_resource"],
        intervals=config["target"],
        gatk=config["gatk"]
    log:
        "logs/mutect2_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} Mutect2 "
        "-R {params.reference} "
        "-I {input.tumor} "
        "-L {params.intervals} "
        "--germline-resource {params.germline} "
        "--f1r2-tar-gz {output.tar} "
        "-O {output.vcf}) 2> {log}"

rule learn_read_orientation_model:
    input:
        "results_tumorOnly/{tumors}/unfiltered.f1r2.tar.gz"
    output:
        protected("results_tumorOnly/{tumors}/read_orientation_model.tar.gz")
    params:
        gatk=config["gatk"]
    log:
        "logs/learn_read_orientation_model_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} LearnReadOrientationModel "
        "-I {input} "
        "-O {output}) 2> {log}"

rule get_pileup_summaries_tumor_only:
    input:
        lambda wc: config["samples"][wc.tumors]
    output:
        protected("results_tumorOnly/{tumors}/pileup_summaries.table")
    params:
        gatk=config["gatk"],
        reference_genome=config["reference_genome"],
        known=config["known_polymorphic_sites"]
    log:
        "logs/get_pileup_summaries_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} GetPileupSummaries "
        "-I {input} "
        "-V {params.known} "
        "-L {params.known} "
        "-R {params.reference_genome} "
        "-O {output}) 2> {log}"

rule calculate_contamination_tumor_only:
    input:
        "results_tumorOnly/{tumors}/pileup_summaries.table"
    output:
        segments=protected("results_tumorOnly/{tumors}/segments.table"),
        contamination=protected("results_tumorOnly/{tumors}/contamination.table")
    params:
        gatk=config["gatk"]
    log:
        "logs/calculate_contamination_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} CalculateContamination "
        "-I {input} "
        "-tumor-segmentation {output.segments} "
        "-O {output.contamination}) 2> {log}"

rule filter_mutect_calls_tumor_only:
    input:
        vcf="results_tumorOnly/{tumors}/unfiltered.vcf.gz",
        tbi="results_tumorOnly/{tumors}/unfiltered.vcf.gz.tbi",
        segments="results_tumorOnly/{tumors}/segments.table",
        contamination="results_tumorOnly/{tumors}/contamination.table",
        model="results_tumorOnly/{tumors}/read_orientation_model.tar.gz",
        stats="results_tumorOnly/{tumors}/unfiltered.vcf.gz.stats"
    output:
        filtered=protected("results_tumorOnly/{tumors}/filtered_all.vcf.gz"),
        filtering_stats=protected("results_tumorOnly/{tumors}/filtering_stats.tsv")
    params:
        gatk=config["gatk"],
        reference=config["reference_genome"]
    log:
        "logs/filter_mutect_calls_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} FilterMutectCalls "
        "-R {params.reference} "
        "-V {input.vcf} "
        "--tumor-segmentation {input.segments} "
        "--contamination-table {input.contamination} "
        "--ob-priors {input.model} "
        "--stats {input.stats} "
        "--filtering-stats {output.filtering_stats} "
        "-O {output.filtered}) 2> {log}"

rule extract_pass_variants:
    input:
        "results_tumorOnly/{tumors}/filtered_all.vcf.gz"
    output:
        protected("results_tumorOnly/{tumors}/pass_variants.vcf.gz")
    params:
        gatk=config["gatk"],
        reference=config["reference_genome"]
    log:
        "logs/extract_pass_variants_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} SelectVariants "
        "-R {params.reference} "
        "-V {input} "
        "--exclude-filtered "
        "--create-output-variant-index "
        "-O {output}) 2> {log}"

rule extract_pass_snvs:
    input:
        "results_tumorOnly/{tumors}/filtered_all.vcf.gz"
    output:
        protected("results_tumorOnly/{tumors}/pass_snvs.vcf.gz")
    params:
        gatk=config["gatk"],
        reference=config["reference_genome"]
    log:
        "logs/extract_pass_snvs_tumor_only/{tumors}.txt"
    shell:
        "({params.gatk} SelectVariants "
        "-R {params.reference} "
        "-V {input} "
        "--select-type-to-include SNP "
        "--exclude-filtered "
        "--create-output-variant-index "
        "-O {output}) 2> {log}"

rule runAnnovar_variants:
    input:
        "results_tumorOnly/{tumors}/pass_variants.vcf.gz"
    output:
        protected("results_tumorOnly/{tumors}/pass_variants.hg38_multianno.vcf")
    params:
        script=config["annovar_python_script"],
        interpreter=config["interpreter"]
    log:
        "logs/runAnnovar_tumor_only/{tumors}_variants.txt"
    shell:
        "({params.interpreter} {params.script} --input_vcf_file_path {input}) 2> {log}"

rule runAnnovar_snvs:
    input:
        "results_tumorOnly/{tumors}/pass_snvs.vcf.gz"
    output:
        protected("results_tumorOnly/{tumors}/pass_snvs.hg38_multianno.vcf")
    params:
        script=config["annovar_python_script"],
        interpreter=config["interpreter"]
    log:
        "logs/runAnnovar_tumor_only/{tumors}_snvs.txt"
    shell:
        "({params.interpreter} {params.script} --input_vcf_file_path {input}) 2> {log}"

rule runAnnovar_filtered:
    input:
        "results_tumorOnly/{tumors}/filtered_all.vcf.gz"
    output:
        protected("results_tumorOnly/{tumors}/filtered_all.hg38_multianno.vcf")
    params:
        script=config["annovar_python_script"],
        interpreter=config["interpreter"]
    log:
        "logs/runAnnovar_tumor_only/{tumors}_filtered.txt"
    shell:
        "({params.interpreter} {params.script} --input_vcf_file_path {input}) 2> {log}"
