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
    --file-prefix {wildcards.local}_{wildcards.window_size}_hudson_ \
    --allow-file-overwriting"  ##pops=$(cat {input.rename} | awk '$2 ~ \"{wildcards.local}\"' | cut -f 2 | sed ':a;N;$!ba;s/\\n/,/g' )
#    --window-average-policy valid-snps \
#-------------------------------------------------------------------------------------------------------------------------------------------------
#locals = list(config['local_to_pairs'].keys())
#all_cutoff_files = [
#    os.path.join(config['analysis_path'], f"fst/window/cutoff_{local}_{ws}_{pair}_{chrom}_.tsv")
#    for local in locals
#    for pair in config['local_to_pairs'][local]
#    for ws in ["100", "250"]
#    for chrom in config['chrom']
#]

rule w_fst_cutoffs:
    input: os.path.join(config['analysis_path'], "fst/window/{local}_{window_size}_hudson_fst.csv"),
    output:
        #manh_wfst = "../results/manh_wfst_{local}_{window_size}.jpeg",
        w_size = os.path.join(config['analysis_path'], "fst/window/BP_windown_size_{local}_{window_size}_{chrom}.tsv"),
        #w_fst_cutoff = lambda wildcards: [
        #    os.path.join(
        #        config['analysis_path'],
        #        f"fst/window/cutoff_{wildcards.local}_{wildcards.window_size}_{pair}_{chrom}_.tsv"
        #    )
        #    for pair in config['local_to_pairs'][wildcards.local]
        #    for chrom in config['chrom']
        #]
    params: w_fst_cutoff = os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)",  window_size = "(100|250)"
    shell: "Rscript scripts/R/outliers_window_FST_per_pop_comp.R -wfst {input} -tfst {params.w_fst_cutoff} -wsize {output.w_size} -chrom {wildcards.chrom}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule merge_top_windows:
    input:
         w_size = os.path.join(config['analysis_path'], "fst/window/BP_windown_size_{local}_{window_size}_{chrom}.tsv"),
         w_fst_cutoff = os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}_{pair}_{chrom}_.tsv")
         
    output: temp(os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_{chrom}.tsv"))
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "awk '$6 <= {wildcards.cutoff}' {input.w_fst_cutoff} | bedtools merge -d 6000 -c 1 -o count > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule get_genes_in_wind:
    input:
        windows = os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_{chrom}.tsv"),
        genes = '../resources/annot_{chrom}.txt'
    output: temp(os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_{chrom}_with_genes.tsv"))
    shell: "python scripts/python/FBgn_assign.py {input.windows} {input.genes} {wildcards.chrom} {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule join_genes_chroms:
    input: expand(os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{{cutoff}}_{{local}}_{{window_size}}_{{pair}}_{chrom}_with_genes.tsv"), chrom = ["2L","2R","3L","3R","X"])
    output: os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes.tsv")
    shell: "cat {input} | awk -v c=\"chrm\" 'NR==1 || index($0,c)==0' > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule enrichment_GO_outliers_fst_windows:
    input:
        w = expand(os.path.join(config['analysis_path'], "fst/window/cutoff_{{local}}_{{window_size}}_{{pair}}_{chrom}_.tsv"), chrom = ["2L","2R","3L","3R","X"]),
        f = os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes.tsv"),
        a = expand('../resources/annot_{chrom}.txt', chrom = ["2L","2R","3L","3R","X"])
    params: 
        w_prefix = os.path.join(config['analysis_path'], "fst/window/cutoff_{local}_{window_size}_{pair}_"),
        a_prefix = '../resources/annot_'
    output: os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes_GO{batch}.txt")
    wildcard_constraints: batch = "|".join([str(i) for i in range(1, 21)])
    shell: "perl scripts/perl/GO_overlap_parallel.pl -f {input.f} -w {params.w_prefix} -s 3 -a {params.a_prefix} -r 500 -b {wildcards.batch}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule merge_GOs_output:
    input: expand(os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{{cutoff}}_{{local}}_{{window_size}}_{{pair}}_all_chroms_with_genes_GO{batch}.txt"), batch = [str(i) for i in range(1, 21)])
    output: os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes_GO_merged.txt")
    params: go_prefix = os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes_GO")
    shell: 'perl scripts/perl/GO_overlap_parallel_merge.pl -i {params.go_prefix}'

#-------------------------------------------------------------------------------------------------------------------------------------------------
#does the same thin as merge_top_windows, but keep some stats
#not the most efficient thing, but all this problem was caused by the perl scripts
rule merge_top_windows_with_stats:
    input:
         w_size = expand(os.path.join(config['analysis_path'], "fst/window/BP_windown_size_{{local}}_{{window_size}}_{chrom}.tsv"), chrom = config['chrom']),
         w_fst_cutoff = expand(os.path.join(config['analysis_path'], "fst/window/cutoff_{{local}}_{{window_size}}_{{pair}}_{chrom}_.tsv"), chrom = config['chrom'])
         
    output: os.path.join(config['analysis_path'], "fst/window/top_{cutoff}_{local}_{window_size}_{pair}_all.tsv")
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "cat {input.w_fst_cutoff} | awk -v c=\"chrm\" 'NR==1 || index($0,c)==0' | awk '$6 <= {wildcards.cutoff}' | bedtools merge -d 6000 -c 1,4,4,5 -o count,max,mean,distinct > {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule join_chrom_and_genes:
    input: 
        top = os.path.join(config['analysis_path'], "fst/window/top_{cutoff}_{local}_{window_size}_{pair}_all.tsv"),
        gft = os.path.join(config['ref_folder'], 'dmel-6.55.gtf')
    output: os.path.join(config['analysis_path'], "fst/window/top_{cutoff}_{local}_{window_size}_{pair}_all_with_gft.tsv")
    shell: "bedtools intersect -wa -wb  -a {input.top} -b {input.gft} | awk '$10 == \"gene\"' > {output}" 

#-------------------------------------------------------------------------------------------------------------------------------------------------
#finds intersect between gft and top windows (already merged together by bedtools merge)
rule intersect_gft:
    input:
        merged_wind = os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes.tsv"),
        gft = os.path.join(config['ref_folder'], 'dmel-6.55.gtf')
    output: os.path.join(config['analysis_path'], "fst/window/merged_cutoffs_{cutoff}_{local}_{window_size}_{pair}_all_chroms_with_genes_gft.tsv")
    shell: "tail -n +2 {input.merged_wind} | cut -f 1,2,3,4 | tee ../resources/temp_file.txt | bedtools intersect -wa -wb  -a - -b {input.gft} > {output}"

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
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)"
    shell: 
        "pair=($(tail -n +2 {input} | sort -k 4,4 | cut -f 4 | uniq )) && "
        "for i in \"${{pair[@]}}\"; do "
            "echo $i ;"
            "awk '$7 == 0.001 && $4 == \"'\"$i\"'\" {{print $0}}' {input} | cut -f 1-5 | sort -k 1,1 -k 2,2n | bedtools merge -d 4000 -c 1,5,5,4 -o count,max,mean,distinct -i - >> {output}; "
        "done"
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule separate_fst_pairs:
    input: expand(os.path.join(config['analysis_path'], "fst/window/top_0.001_{local}_100.bed"), local = ["MA","CT","VT","MD","MFL","HFL"]),
    output: os.path.join(config['analysis_path'], "fst/window/GOwinda/candidate_windows_{pair}.tsv")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)", pair = "(CMA97\\.HMA09|MCT97\\.MCT09|dlFL10\\.HFL97|SVT09\\.WVT97|CMD10\\.CMD97B|CMA97\\.MA22|MCT97\\.CT22|MFL97\\.MFL23)"
    shell: "grep {wildcards.pair} {input} | cut -d \":\" -f 2 > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule get_snps_within_outliers:
    input: 
        top_w = os.path.join(config['analysis_path'], "fst/window/GOwinda/candidate_windows_{pair}.tsv"),
        vcf = os.path.join(config['call_path'],"PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: os.path.join(config['analysis_path'], "fst/window/GOwinda/candidate_SNPs_{pair}.tsv")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)", pair = "(CMA97\\.HMA09|MCT97\\.MCT09|dlFL10\\.HFL97|SVT09\\.WVT97|CMD10\\.CMD97B|CMA97\\.MA22|MCT97\\.CT22|MFL97\\.MFL23)"
    shell: "bedtools intersect -wb -a {input.top_w} -b {input.vcf} | cut -f 8,9 > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule GO_wfst:
    input:
        gtf=os.path.join(config['ref_folder'], "dmel-6.55.gtf"),
        goassociation =  os.path.join(config['ref_folder'],"GO_gene_association.tsv"),
        total_snps = os.path.join(config['call_path'],"PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        cand_snps= os.path.join(config['analysis_path'], "fst/window/GOwinda/candidate_SNPs_{pair}.tsv")
    output: os.path.join(config['analysis_path'], "fst/window/GOwinda/GO_enrichment_mode_gene_gendef_{gen_def}_pair_{pair}_top_0.001_fst_windows.tsv")
    threads: 8
    params: output_name="GO_enrichment_mode_gene_gendef_{gen_def}_pair_{pair}_top_0.001_fst_windows.tsv"
    wildcard_constraints: gen_def = "updownstream1000", pair = "(CMA97\\.HMA09|MCT97\\.MCT09|dlFL10\\.HFL97|SVT09\\.WVT97|CMD10\\.CMD97B|CMA97\\.MA22|MCT97\\.CT22|MFL97\\.MFL23)"
    shell: "java -Xmx4g -jar /home/vitoria/bin/Gowinda-1.12.jar --snp-file {input.total_snps} \
    --candidate-snp-file {input.cand_snps} \
    --gene-set-file {input.goassociation} \
    --annotation-file {input.gtf} \
    --simulations 100000 \
    --min-significance 1 \
    --gene-definition {wildcards.gen_def} \
    --threads {threads} \
    --output-file {params.output_name} \
    --mode gene \
    --min-genes 5 || mv {params.output_name} {output}"


#-------------------------------------------------------------------------------------------------------------------------------------------------
rule merge_outliers_windows_with_genes: 
    input: 
        top_w = os.path.join(config['analysis_path'], "fst/window/top_0.001_{local}_{window_size}.bed"),
        gtf = os.path.join(config['ref_folder'], "dmel-6.55.gtf")
    output: os.path.join(config['analysis_path'], "fst/window/top_0.001_{local}_{window_size}_with_gtf.bed")
    wildcard_constraints: local = "(MA|CT|VT|MD|MFL|HFL)"
    shell: "bedtools intersect -wa -wb -header -a {input.top_w} -b {input.gtf} | awk '$10 == \"gene\"'> {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule bam_txt_cyp:
    input: 
        bams = expand(os.path.join(config['align_path'], '{cyp_pops}.md.srt.bam'), cyp_pops = ["HFL97downto60mi_L001", "dlFL10_L001", "MFL23_L002","12_L001", "CMA97_L002","HMA09_L002","A45_L002", "MCT09_L002", "MCT97_L002", "A46_L002"]), 
        bais = expand(os.path.join(config['align_path'], '{cyp_pops}.md.srt.bam.bai'), cyp_pops = ["HFL97downto60mi_L001", "dlFL10_L001", "MFL23_L002","12_L001", "CMA97_L002","HMA09_L002","A45_L002", "MCT09_L002", "MCT97_L002", "A46_L002"]),
    output: os.path.join(config['align_path'], "bam_list_cyp.txt")
    #wildcard_constraints: id="[a-zA-Z0-9]+", lane="\d"
    shell:
        "echo {input.bams} | sed 's/ /\\n/g' > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# - computes depth per position
rule compute_depth_stats_cyp:
    input: 
        bam_list_txt= os.path.join(config['align_path'], "bam_list_cyp.txt"), 
        bams = expand(os.path.join(config['align_path'], '{cyp_pops}.md.srt.bam'), cyp_pops = ["HFL97downto60mi_L001", "dlFL10_L001", "MFL23_L002","12_L001", "CMA97_L002","HMA09_L002","A45_L002", "MCT09_L002", "MCT97_L002", "A46_L002"]),
        bais = expand(os.path.join(config['align_path'], '{cyp_pops}.md.srt.bam.bai'), cyp_pops = ["HFL97downto60mi_L001", "dlFL10_L001", "MFL23_L002","12_L001", "CMA97_L002","HMA09_L002","A45_L002", "MCT09_L002", "MCT97_L002", "A46_L002"]),
    output: os.path.join(config['qltctrl_path'], "pre_downsample/samtools_depth_stats_cyp_depth_include_unmerged_flag.tsv")
    #wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell: "samtools depth -r 2R:14800000-14900000 -H -f {input.bam_list_txt} -q 20 -Q 20 -o {output}"
#cyp_depths_plot.R para gerar o plot
#cyp_depths_plot.pdf para ver a figura
#-------------------------------------------------------------------------------------------------------------------------------------
#próxima regra é um Rscript que pega as janelas sobrepostas entre diferentes localidades para uma mesma comparação temporal (o pré-script está na pasta clinas "shared_fst_outliers_space.R")
#seguida de bedtools merge que pega as genes que se sobrepõem à janela do gtf
#bedtools merge -d 10000 -c 5,6,7,8,13 -o distinct -i shared_windows_with_genes.tsv
#por último virá a regra que vai gerar a figura com as regiões genomicas com janelas sobrepostas e os arquivos com os genes e gos (o pré-script está na pasta clinas "shered_genes_fst_outliers.R")