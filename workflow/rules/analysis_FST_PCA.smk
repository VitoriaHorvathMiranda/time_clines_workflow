
##PRECISA MUDAR ESSA FUNÇÃO PARA UM CHECKPOINT NO FUTURO
def get_pairs_from_file(file):
    pairs = set()
    with open(file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 5:  # Ensure there are at least 5 columns
                pairs.add(columns[4])
    pairs = list(pairs)
    pairs.sort()  # Sort the pairs if needed
    return pairs
#-------------------------------------------------------------------------------------------------------------------------------------------------

rule PCA:
    input:
        meta=config['meta_path'],
        freqs=os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    output:
        PC1x2all= "../results/PC1x2_all_pops.jpeg",
        PC3x4all= "../results/PC3x4_all_pops.jpeg",
        PCAtime= "../results/PCA_per_time.jpeg",
        PCAchrom="../results/PCA_per_chrom.jpeg",
        lmALL = os.path.join(config['analysis_path'], "PCA/lm_PC1-4_all_variables.tsv"),
        Rsq = os.path.join(config['analysis_path'], "PCA/Rsq_PC1-4_all_variables.tsv")
    shell: "Rscript scripts/R/PCA_script.R -m {input.meta} -f {input.freqs} -p12All {output.PC1x2all} -p34All {output.PC3x4all} -pTimes {output.PCAtime} -PCAchrom {output.PCAchrom} -lm {output.lmALL} -Rsq {output.Rsq}"

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
# filter problematic regions out of windowns to get true window size;
# martin Kapun Script from https://github.com/capoony/DrosEU_pipeline/tree/master #small uptate to accept .bed file instead of .gff 
rule true_windowns:
    input:
        bad_pos = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_BS.txt.gz"),
        indels = os.path.join(config['align_path'], "InDel-positions_20.txt.gz"),
        te = os.path.join(config['ref_folder'], "repeat_6_with_spaces.bed")
    output: os.path.join(config['analysis_path'], "pop_stats/truewindows-200000-200000.txt"), 
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/truewindows-200000-200000")
    shell: "python2 scripts/python/TrueWindows.py \
    --badcov {input.bad_pos} \
    --indel {input.indels} \
    --te {input.te} \
    --window 200000 \
    --step 200000 \
    --chromosomes 2L:23513712,2R:25286936,3L:28110227,3R:32079331  \
    --output {params.out_path}" #,X:23542271

#-------------------------------------------------------------------------------------------------------------------------------------------------
#computes thetaPi, thetaW and tajima's D 
# martin Kapun Script from https://github.com/capoony/DrosEU_pipeline/tree/master
rule pop_stats:
    input: 
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows-200000-200000.txt"),
        pool_sizes = "../resources/pool_sizes_autosome.csv"
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.th")
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000")
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --pool-size $pool_sizes \
    --min-sites-frac 0.6 \
    --window 200000 \
    --step 200000 \
    --sitecount {input.sitecount} \
    --output {params.out_path}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
# filter problematic regions out of windowns to get true window size;
# martin Kapun Script from https://github.com/capoony/DrosEU_pipeline/tree/master #small uptate to accept .bed file instead of .gff 
rule true_windowns_X:
    input:
        bad_pos = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_BS.txt.gz"),
        indels = os.path.join(config['align_path'], "InDel-positions_20.txt.gz"),
        te = os.path.join(config['ref_folder'], "repeat_6_with_spaces.bed")
    output: os.path.join(config['analysis_path'], "pop_stats/truewindows_X-200000-200000.txt"), 
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/truewindows_X-200000-200000")
    shell: "python2 scripts/python/TrueWindows.py \
    --badcov {input.bad_pos} \
    --indel {input.indels} \
    --te {input.te} \
    --window 200000 \
    --step 200000 \
    --chromosomes X:23542271  \
    --output {params.out_path}" #,

#-------------------------------------------------------------------------------------------------------------------------------------------------
#computes thetaPi, thetaW and tajima's D 
# martin Kapun Script from https://github.com/capoony/DrosEU_pipeline/tree/master
rule pop_stats_X:
    input: 
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows_X-200000-200000.txt"),
        pool_sizes = "../resources/pool_sizes_X.csv"
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.th")
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000")
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --pool-size $pool_sizes \
    --min-sites-frac 0.6 \
    --window 200000 \
    --step 200000 \
    --sitecount {input.sitecount} \
    --output {params.out_path}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 
rule plot_stats:
    input:
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf"),
        meta = config['meta_path'],
        a = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.pi"), b = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.D"), c = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.th"),
        e = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.pi"), f = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.D"), g = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.th")
    output: "../results/pop_stats_pi.jpeg", "../results/pop_stats_w.jpeg", "../results/pop_stats_D.jpeg"
    params: stats_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000") , plot_path = "../results/pop_stats_"
    shell: "Rscript scripts/R/pop_stats_kapun.R --vcf {input.vcf} --meta {input.meta} --stats {params.stats_path} --plots {params.plot_path}"


#-------------------------------------------------------------------------------------------------------------------------------------------------
rule pop_stats_chrom_arms:
    input:
        statA = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_miss_fraq0.60_200000_200000.{stat}"),
        statX = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_miss_fraq0.60_200000_200000.{stat}"),
        windowA = os.path.join(config['analysis_path'], "pop_stats/truewindows-200000-200000.txt"),
        windowX = os.path.join(config['analysis_path'], "pop_stats/truewindows_X-200000-200000.txt"),
        meta = config['meta_path'],
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output:
        per_pop = os.path.join(config['analysis_path'], "pop_stats/Global_stat_{stat}_per_pop.tsv"),
        per_year = os.path.join(config['analysis_path'], "pop_stats/Global_stat_{stat}_per_year.tsv"),
        barplot = "../results/barplot_stat_{stat}.jpeg"
    wildcard_constraints: stat = "(pi|th)"
    shell: "Rscript scripts/R/global_pop_stat.R -statA {input.statA} -statX {input.statX} \
    -wA {input.windowA} -wX {input.windowX} \
    -meta {input.meta} -vcf {input.vcf} \
    -out {output.per_pop} -y {output.per_year} -bar {output.barplot}"

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
rule FST_per_chrom:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes_autosome.csv",
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv"
    output:fst = os.path.join(config['FST_genome_path'],"Genome_FST_{chrom}_fst-list.csv"),
    params: outdir=os.path.join(config['analysis_path'], "fst/genome"), prefix="Genome_FST_{chrom}_"
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "/home/vitoria/bin/grenedalfv6.0/grenedalf-0.6.0/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta {input.ref} \
    --filter-region {wildcards.chrom}\
    --window-type genome \
    --method unbiased-hudson \
    --rename-samples-list {input.rename} \
    --filter-sample-min-read-depth 10 \
    --pool-sizes {input.pool_sizes} \
    --window-average-policy valid-snps \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule X_FST:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        pool_sizes = "../resources/pool_sizes_X.csv",
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv"
    output: os.path.join(config['FST_genome_path'],"Genome_FST_X_fst-matrix.csv"), os.path.join(config['FST_genome_path'],"Genome_FST_X_fst-list.csv")
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
rule fst_per_year:
    input: 
        fst = os.path.join(config['FST_genome_path'],"Genome_FST_{chrom_type}_fst-list.csv"),
        meta = config['meta_path']
    output:
        boxplot = "../results/fst_year_boxplot_{chrom_type}.jpeg",
        lm = "../results/fst_year_lm_{chrom_type}.txt",
        distLM = "../results/fst_year_dist_lm_{chrom_type}.txt",
        distfig = "../results/fst_year_dist_{chrom_type}.jpeg",
    wildcard_constraints: chrom_type = "(autosome|X|2L|2R|3L|3R)"
    shell: "Rscript scripts/R/mean_fst_per_year.R -fst {input.fst} -meta {input.meta} -lm {output.lm} -box {output.boxplot} -distLM {output.distLM} -dist {output.distfig}"

#-------------------------------------------------------------------------------------------------------------------------------------------------

#rule


#rule thetaPI:
#    input:
#        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
#        pool_sizes = "../resources/pool_sizes_autosome.csv",
#        rename="../resources/rename_samples.tsv"
#    output: temp(os.path.join(config['thetaPI_path'], "thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv"))
#    params:  outdir=config['thetaPI_path'], prefix="thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97"
#    threads: 40
#    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf diversity \
#    --sync-path {input.sync} \
#    --filter-sample-min-count 1 \
#    --filter-sample-min-coverage 15 \
#    --window-type single \
#    --measure theta-pi \
#    --rename-samples-file {input.rename} \
#    --pool-sizes {input.pool_sizes} \
#    --out-dir {params.outdir} \
#    --file-prefix {params.prefix} \
#    --threads {threads}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#prestar atenção no número de pops e no número de grupos '0,0.000,0,0'
#rule filter_thetaPI:
#    input: os.path.join(config['thetaPI_path'], "thetaPi_single_with_dlGA10_dlSC10_noSNC10_noESC97diversity.csv")
#    output: os.path.join(config['thetaPI_path'], "thetaPi_filtered.csv")
#    shell: "grep -v '0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0,0,0.000,0,0' {input} > {output}"


#rule true_windowns_chrom_arm:
#    input:
#        bad_pos = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_BS.txt.gz"),
#        indels = os.path.join(config['align_path'], "InDel-positions_20.txt.gz"),
#        te = os.path.join(config['ref_folder'], "repeat_6_with_spaces.bed")
#    output: os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
#    params: 
#        out_path = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}"),
#        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1],
#        chrom_value = lambda wildcards: f"{wildcards.chrom_lenght_step.split('-')[0]}:{wildcards.chrom_lenght_step.split('-')[1]}"
#    shell: "python2 scripts/python/TrueWindows.py \
#    --badcov {input.bad_pos} \
#    --indel {input.indels} \
#    --te {input.te} \
#    --window {params.window_size} \
#    --step {params.window_size} \
#    --chromosomes {params.chrom_value} \
#    --output {params.out_path}" 

#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule pop_stats_global:
#    input:
#        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
#        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
#        pool_sizes = "../resources/pool_sizes_autosome.csv"
#    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.th")
#    wildcard_constraints: chrom_lenght_step = "|".join([v for v in config['chrom_lenght_step'] if v != "X-23542271-23542271"])
#    params: 
#        out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}"),
#        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1]
#    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
#    python2 scripts/python/PoolGen-var.py \
#    --input {input.sync} \
#    --min-count 1 \
#    --min-sites-frac 0.5 \
#    --pool-size $pool_sizes \
#    --window {params.window_size} \
#    --step {params.window_size} \
#    --sitecount {input.sitecount} \
#    --output {params.out_path}"
        
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule pop_stats_global_X:
#    input:
#        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
#        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
#        pool_sizes = "../resources/pool_sizes_X.csv"
#    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_global_{chrom_lenght_step}.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_global_{chrom_lenght_step}.th")
#    wildcard_constraints: chrom_lenght_step = ("X-23542271-23542271")
#    params: 
#        out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}"),
#        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1]
#    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
#    python2 scripts/python/PoolGen-var.py \
#    --input {input.sync} \
#    --min-count 1 \
#    --min-sites-frac 0.5 \
#    --pool-size $pool_sizes \
#    --window {params.window_size} \
#    --step {params.window_size} \
#    --sitecount {input.sitecount} \
#    --output {params.out_path}"

