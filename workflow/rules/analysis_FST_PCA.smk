


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
    input: 
        meta = config['meta_path'], 
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: 
        auto = "../resources/pool_sizes_autosome.csv", 
        X = "../resources/pool_sizes_X.csv",
        rename="../resources/rename_samples.tsv"
    params: sync_name="PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h"
    shell: "Rscript scripts/R/make_pool_sizes_file.R -meta {input.meta} -vcf {input.vcf} -string {params.sync_name} -oa {output.auto} -ox {output.X} -on {output.rename}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule thetaPI:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes_autosome.csv",
        rename="../resources/rename_samples.tsv"
    output: temp(os.path.join(config['thetaPI_path'], "thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv"))
    params:  outdir=config['thetaPI_path'], prefix="thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97"
    threads: 40
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf diversity \
    --sync-path {input.sync} \
    --filter-sample-min-count 1 \
    --filter-sample-min-coverage 15 \
    --window-type single \
    --measure theta-pi \
    --rename-samples-file {input.rename} \
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
        #vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: "../results/thetapi_plot.jpeg"
    shell: "Rscript scripts/R/thetaPI_stat.R -pi {input.pi} -plot {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule global_FST_autosome:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes_autosome.csv",
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv"
    output: os.path.join(config['FST_genome_path'],"Genome_FST_autosome_fst-matrix.csv")
    params: outdir=config['FST_genome_path'], prefix="Genome_FST_autosome_"
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --filter-region 2L --filter-region 2R --filter-region 3L --filter-region 3R \
    --window-type genome \
    --method unbiased-hudson \
    --rename-samples-file {input.rename} \
    --pool-sizes {input.pool_sizes} \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule X_FST:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes_X.csv",
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv"
    output: os.path.join(config['FST_genome_path'],"Genome_FST_X_fst-matrix.csv")
    params: outdir=config['FST_genome_path'], prefix="Genome_FST_X_"
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --filter-region X \
    --window-type genome \
    --method unbiased-hudson \
    --rename-samples-file {input.rename} \
    --pool-sizes {input.pool_sizes} \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule plot_genome_FST:
    input:
        meta = config['meta_path'],
        fst_x = os.path.join(config['FST_genome_path'],"Genome_FST_X_fst-matrix.csv"),
        fst_auto = os.path.join(config['FST_genome_path'],"Genome_FST_autosome_fst-matrix.csv")
    output:
        all_pops = "../results/ALL_pops_genome_FST.jpeg",
        time_pops = "../results/Time_pops_genome_FST.jpeg"
    shell: "Rscript scripts/R/genome_FST_plot.R -meta {input.meta} -x {input.fst_x} -auto {input.fst_auto} -oall {output.all_pops} -otime {output.time_pops}"



#-------------------------------------------------------------------------------------------------------------------------------------------------
rule window_fst:
    input: 
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv",
        regions = "../resources/regions.bed", #*#regions with recombination above 0.5 cM/Mb - from pool et al 2017 - lifted over for ref 6 genome with https://flybase.org/convert/coordinates
        pool_sizes = "../resources/pool_sizes_autosome.csv" #script que fez esse arquivo está na pasta grenedalf_analysis
    output: os.path.join(config['analysis_path'], "fst/window/{local}_{window_size}_fst.csv")
    params: out_dir = os.path.join(config['analysis_path'], "fst/window/")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL)", window_size = "(100|250)"
    shell: "pops=$(cat {input.rename} | awk '$2 ~ \"{wildcards.local}\"' | cut -f 2 | sed ':a;N;$!ba;s/\\n/,/g' ) && \
    stride=$(awk 'BEGIN {{ printf \"%.0f\\n\", 0.1 * {wildcards.window_size}}}') && \
    /home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --rename-samples-file {input.rename} \
    --filter-samples-include $pops \
    --filter-region-bed {input.regions} \
    --window-type queue \
    --window-queue-count {wildcards.window_size} \
    --window-queue-stride $stride \
    --method unbiased-nei \
    --pool-sizes {input.pool_sizes} \
    --verbose \
    --out-dir {params.out_dir} \
    --file-prefix {wildcards.local}_{wildcards.window_size}_"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule SNP_fst:
    input: 
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv",
        regions = "../resources/regions.bed", #*#regions with recombination above 0.5 cM/Mb - from pool et al 2017 - lifted over for ref 6 genome with https://flybase.org/convert/coordinates
        pool_sizes = "../resources/pool_sizes_autosome.csv" #script que fez esse arquivo está na pasta grenedalf_analysis
    output: os.path.join(config['analysis_path'], "fst/SNPs/{local}_SNPs_fst.csv")
    params: out_dir = os.path.join(config['analysis_path'], "fst/SNPs/")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL)"
    shell: "pops=$(cat {input.rename} | awk '$2 ~ \"{wildcards.local}\"' | cut -f 2 | sed ':a;N;$!ba;s/\\n/,/g' ) && \
    /home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --rename-samples-file {input.rename} \
    --filter-samples-include $pops \
    --filter-region-bed {input.regions} \
    --window-type single \
    --omit-na-windows \
    --method unbiased-nei \
    --pool-sizes {input.pool_sizes} \
    --verbose \
    --out-dir {params.out_dir} \
    --file-prefix {wildcards.local}_SNPs_"



