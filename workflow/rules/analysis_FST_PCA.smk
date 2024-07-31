
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
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.th")
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000")
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --pool-size $pool_sizes \
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
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.th")
    params: out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000")
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --pool-size $pool_sizes \
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
        a = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.pi"), b = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.D"), c = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_X_200000_200000.th"),
        e = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.pi"), f = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.D"), g = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000.th")
    output: "../results/pop_stats_pi.jpeg", "../results/pop_stats_w.jpeg", "../results/pop_stats_D.jpeg"
    params: stats_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_auto_200000_200000") , plot_path = "../results/pop_stats_"
    shell: "Rscript scripts/R/pop_stats_kapun.R --vcf {input.vcf} --meta {input.meta} --stats {params.stats_path} --plots {params.plot_path}"


#-------------------------------------------------------------------------------------------------------------------------------------------------

rule true_windowns_chrom_arm:
    input:
        bad_pos = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_BS.txt.gz"),
        indels = os.path.join(config['align_path'], "InDel-positions_20.txt.gz"),
        te = os.path.join(config['ref_folder'], "repeat_6_with_spaces.bed")
    output: os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
    params: 
        out_path = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}"),
        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1],
        chrom_value = lambda wildcards: f"{wildcards.chrom_lenght_step.split('-')[0]}:{wildcards.chrom_lenght_step.split('-')[1]}"
    shell: "python2 scripts/python/TrueWindows.py \
    --badcov {input.bad_pos} \
    --indel {input.indels} \
    --te {input.te} \
    --window {params.window_size} \
    --step {params.window_size} \
    --chromosomes {params.chrom_value} \
    --output {params.out_path}" 

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule pop_stats_global:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
        pool_sizes = "../resources/pool_sizes_autosome.csv"
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.th")
    wildcard_constraints: chrom_lenght_step = "|".join([v for v in config['chrom_lenght_step'] if v != "X-23542271-23542271"])
    params: 
        out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}"),
        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1]
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --min-count 1 \
    --min-sites-frac 0.5 \
    --pool-size $pool_sizes \
    --window {params.window_size} \
    --step {params.window_size} \
    --sitecount {input.sitecount} \
    --output {params.out_path}"
        
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule pop_stats_global_X:
    input:
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        sitecount = os.path.join(config['analysis_path'], "pop_stats/truewindows_chrom_arm_{chrom_lenght_step}.txt"),
        pool_sizes = "../resources/pool_sizes_X.csv"
    output: os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}.pi"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_global_{chrom_lenght_step}.D"), os.path.join(config['analysis_path'], "pop_stats/pop_stats_kapun_global_{chrom_lenght_step}.th")
    wildcard_constraints: chrom_lenght_step = ("X-23542271-23542271")
    params: 
        out_path = os.path.join(config['analysis_path'], "pop_stats/pop_stats_global_{chrom_lenght_step}"),
        window_size = lambda wildcards: wildcards.chrom_lenght_step.split('-')[1]
    shell: "pool_sizes=$(cut -d ',' -f 2 {input.pool_sizes} | paste -sd,) && \
    python2 scripts/python/PoolGen-var.py \
    --input {input.sync} \
    --min-count 1 \
    --min-sites-frac 0.5 \
    --pool-size $pool_sizes \
    --window {params.window_size} \
    --step {params.window_size} \
    --sitecount {input.sitecount} \
    --output {params.out_path}"

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
rule window_fst: #não precisa fazer uma regra separada para o cromossomo X porque todas as populações que eu uso aqui são só de fêmeas e foram coletadas da mesma maneira
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
rule w_fst_cutoffs:
    input: os.path.join(config['analysis_path'], "fst/window/{local}_{window_size}_fst.csv"),
    output:
        manh_wfst = "../results/manh_wfst_{local}_{window_size}.jpeg",
        w_size = os.path.join(config['analysis_path'], "fst/window/BP_windown_size_{local}_{window_size}.tsv"),
        w_fst_cutoff = os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}.tsv")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL)", window_size = "(100|250)"
    shell: "Rscript scripts/R/outliers_window_FST_all_comp.R -wfst {input} -tfst {output.w_fst_cutoff} -wsize {output.w_size} -fig {output.manh_wfst}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule SNP_fst: #não precisa fazer uma regra separada para o cromossomo X porque todas as populações que eu uso aqui são só de fêmeas e foram coletadas da mesma maneira
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

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule snp_fst_cutoffs:
    input: os.path.join(config['analysis_path'], "fst/SNPs/{local}_SNPs_fst.csv")
    output: 
        manh_wfst = "../results/manh_snp_{local}.jpeg",
        snp_fst_cutoff = os.path.join(config['analysis_path'], "fst/SNP/cutoff_{local}_snp.tsv")
    shell: "Rscript scripts/R/outliers_snp_FST_all_comp.R -sfst {input} -tfst {output.snp_fst_cutoff} -fig {output.manh_wfst}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#tem que ver se a função do começo desse .smk e o uso dela na rule all vai funcionar
rule separate_fst_pairs:
    input: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}.tsv")
    output: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}_{pair}.tsv")
    shell:"awk '$5 == \"{wildcards.pair}\" {{print $0}}' {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# filtra par aos top 0.1% das janelas e funde as janelas com sobreposição; fica com o número de janelas fusionadas e o máximo fst da maior janela
rule merge_outliers_windows: 
    input: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}_{pair}.tsv")
    output: os.path.join(config['analysis_path'], "fst/window/top_0.001_{local}_{window_size}_{pair}.bed")
    shell: "awk '$8 == 0.001 {{print $0}}' {input} | cut -f 1,2,3,6 | sort -k 1,1 -k2,2n | bedtools merge -d 1000 -c 1,4 -o count,max -i - > {output}"

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



