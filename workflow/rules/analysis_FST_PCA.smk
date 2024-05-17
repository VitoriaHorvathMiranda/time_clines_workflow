


#-------------------------------------------------------------------------------------------------------------------------------------------------

rule PCA:
    input:
        meta=config['meta_path'],
        freqs=os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    output:
        PCAall= "../results/PCA_all_pops.jpeg",
        PCAtime= "../results/PCA_per_time.jpeg",
        PCAchrom97="../results/PCA_per_chrom_97.jpeg",
        PCAchrom0910="../results/PCA_per_chrom_0910.jpeg"
    shell: "Rscript scripts/R/PCA_script.R -m {input.meta} -f {input.freqs} -pAll {output.PCAall} -pTimes {output.PCAtime} -PC97 {output.PCAchrom97} -PC0910 {output.PCAchrom0910}"


#-------------------------------------------------------------------------------------------------------------------------------------------------
rule pool_sizes:
    input: meta = config['meta_path'], vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: "../resources/pool_sizes.csv"
    shell: "Rscript scripts/R/make_pool_sizes_file.R -meta {input.meta} -vcf {input.vcf} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule thetaPI:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes.csv"
    output: temp(os.path.join(config['thetaPI_path'], "thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv"))
    params:  outdir=config['thetaPI_path'], prefix="thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97"
    threads: 40
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf diversity \
    --sync-path {input.sync} \
    --filter-sample-min-count 1 \
    --filter-sample-min-coverage 15 \
    --window-type single \
    --measure theta-pi \
    --pool-sizes {input.pool_sizes} \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix} \
    --threads {threads}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule filter_thetaPI:
    input: os.path.join(config['thetaPI_path'], "thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv")
    output: os.path.join(config['thetaPI_path'], "thetaPi_filtered.csv")
    shell: "grep -v '0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0' {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule plot_thetaPI:
    input:
        pi = os.path.join(config['thetaPI_path'], "thetaPi_filtered.csv"),
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: "../results/thetapi_plot.jpeg"
    shell: "Rscript scripts/R/thetaPI_stat.R -pi {input.pi} -vcf {input.vcf} -plot {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule global_FST:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes.csv"
    output:
    params: 
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --filter-region-bed {input.regions} \
    --window-type genome \
    --method unbiased-unbiased-hudson \
    --pool-sizes {input.pool_sizes} \
    --out-dir "

