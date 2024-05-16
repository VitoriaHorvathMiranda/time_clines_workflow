


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
rule theta_Pi:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes.csv"
    output: os.path.join(config['thetaPI_path'], "TajimaD_thetaW_thetaPi_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv")
    params:  outdir=config['thetaPI_path'], prefix="TajimaD_thetaW_thetaPi_with_dlGA10_dlSC10_noSNC10_noESC97"
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf diversity \
    --sync-path {input.sync} \
    --filter-sample-min-count 1 \
    --filter-sample-min-coverage 15 \
    --window-type sliding \
    --window-sliding-width 1000 \
    --window-sliding-stride 500 \
    --measure theta-pi \
    --pool-sizes {input.pool_sizes} \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}"