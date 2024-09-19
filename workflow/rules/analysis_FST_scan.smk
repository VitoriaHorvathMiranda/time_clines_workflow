#-------------------------------------------------------------------------------------------------------------------------------------------------
rule window_fst: #não precisa fazer uma regra separada para o cromossomo X porque todas as populações que eu uso aqui são só de fêmeas e foram coletadas da mesma maneira
    input: #talvez precise fazer uma regra para o X por causa dos HFLs
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        ref = config['ref_path'],
        rename = "../resources/rename_samples.tsv",
        regions = "../resources/regions.bed", #*#regions with recombination above 0.5 cM/Mb - from pool et al 2017 - lifted over for ref 6 genome with https://flybase.org/convert/coordinates
        pool_sizes = "../resources/pool_sizes_autosome.csv" #script que fez esse arquivo está na pasta grenedalf_analysis
    output: os.path.join(config['analysis_path'], "fst/window/{local}_{window_size}_hudson_fst.csv")
    params: out_dir = os.path.join(config['analysis_path'], "fst/window/")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)", window_size = "(100|250)"
    shell: "pops=$( [ \"{wildcards.local}\" == \"HFL\" ] && echo \"HFL97,dlFL10\" || awk -v str=\"{wildcards.local}\" '$2 ~ str {{print $2}}' {input.rename} | sed ':a;N;$!ba;s/\\n/,/g' ) && \
    stride=$(awk 'BEGIN {{ printf \"%.0f\\n\", 0.1 * {wildcards.window_size}}}') && \
    /home/vitoria/bin/grenedalfv6.0/grenedalf-0.6.0/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta {input.ref} \
    --rename-samples-list {input.rename} \
    --filter-samples-include $pops \
    --filter-region-bed {input.regions} \
    --window-type queue \
    --window-queue-count {wildcards.window_size} \
    --window-queue-stride $stride \
    --window-average-policy valid-snps \
    --method unbiased-hudson \
    --filter-sample-min-read-depth 10 \
    --pool-sizes {input.pool_sizes} \
    --verbose \
    --out-dir {params.out_dir} \
    --file-prefix {wildcards.local}_{wildcards.window_size}_hudson_" #pops=$(cat {input.rename} | awk '$2 ~ \"{wildcards.local}\"' | cut -f 2 | sed ':a;N;$!ba;s/\\n/,/g' )

#    --window-average-policy valid-snps \
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule w_fst_cutoffs:
    input: os.path.join(config['analysis_path'], "fst/window/{local}_{window_size}_hudson_fst.csv"),
    output:
        manh_wfst = "../results/manh_wfst_{local}_{window_size}.jpeg",
        w_size = os.path.join(config['analysis_path'], "fst/window/BP_windown_size_{local}_{window_size}.tsv"),
        w_fst_cutoff = os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}.tsv")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)", window_size = "(100|250)"
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
#rule separate_fst_pairs:
#    input: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}.tsv")
#    output: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}_{pair}.tsv")
#    shell:"pair=($(sort -k 4,4 {input} | unique )) && awk '$5 == \"{wildcards.pair}\" {{print $0}}' {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# filtra par aos top 0.1% das janelas e funde as janelas com sobreposição; fica com o número de janelas fusionadas e o máximo fst da maior janela
rule merge_outliers_windows: 
    input: os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}.tsv")
    output: os.path.join(config['analysis_path'], "fst/window/top_0.001_{local}_{window_size}.bed")
    shell: 
        "pair=($(tail -n +2 {input} | sort -k 4,4 | cut -f 4 | uniq )) && "
        "for i in \"${{pair[@]}}\"; do "
            "echo $i ;"
            "awk '$7 == 0.001 && $4 == \"'\"$i\"'\" {{print $0}}' {input} | cut -f 1-5 | sort -k 1,1 -k 2,2n | bedtools merge -d 1000 -c 1,5,5,4 -o count,max,mean,distinct -i - >> {output}; "
        "done"
