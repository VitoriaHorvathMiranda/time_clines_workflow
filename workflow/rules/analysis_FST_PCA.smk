


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