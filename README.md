# mutect2_snakemake
A snakemake to run Mutect2 on analysis-ready bams, following GATK best practices. 

This repository contains all the folders needed to run the snakemake. config/ contains configuration files, while logs/ and results/ just contain placeholder files so that the file structure needed to run the snakemake already exists.

To run the snakemake, update config/samples.yaml and config/config.yaml as per the directions in those files, and then follow the instructions in mutect2.snakefile. The snakemake steps will be launched as jobs to the cluster at the Fred Hutchinson Cancer Research Center.

Note: This snakemake can easily be modified to run a tumor-only analysis.  It can also easily be modified to expect NCBI chromosome naming (i.e. "1") instead of UCSC ("chr1"), which it currently expects.  Please contact Anna if you are interested in either of these adjustments until she adds them to the snakemake. :)
