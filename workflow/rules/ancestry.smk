


#-------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------
# the pipeline of the mapping and SNP calling for this vcf is available in the new_samples folder, as is the table detailing the populations used
# 23 central-west African and 22 European (French and Swedish) haploid samples were used 
rule ancestral_painel_prep:
    input: os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.biallelic.clean.vcf"),
    output: os.path.join(config['analysis_path'], "ancestry/raw_ancestry_painel_chrom_{chrom}.tsv")
    shell: "Rscript scripts/R/raw_ancestral_painel_prep.R -vcf {input} -out {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------
rule make_ancestral_sync:
    input: os.path.join(config['analysis_path'], "ancestry/raw_ancestry_painel_chrom_{chrom}.tsv"),
    output: os.path.join(config['analysis_path'], "ancestry/ancestral_chrom_{chrom}.sync")
    shell: "Rscript scripts/R/make_ancestral_sync.R --painel {input} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# computes global ancestry based on bergland et al 2016 paper
# ref ancestral painel snps were filtered for a minimum number of genomes in each painel (European and African) of 15;
# minimal freq difference between ancestral pops of 0.2 (this guarantees the we are looking at informative SNPs only - the same was done in the ahmm paper)
rule global_ancestry:
    input: 
        freq = os.path.join(config['call_path'], "NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv"),
        meta = config['meta_path'],
        painel = os.path.join(config['analysis_path'], "ancestry/raw_ancestry_painel.tsv")
    output:
        allmod = os.path.join(config['analysis_path'], "ancestry/global/all_models_coeff.tsv"),
        figAll = "../results/all_global_ancestry_models_violin.jpeg",
        summary = os.path.join(config['analysis_path'], "ancestry/global/summary_models.tsv"),
        figLat = "../results/global_ancestry_per_lat.jpeg"
    shell: "Rscript scripts/R/ancestry_lm.R -ne {input.freq} \
    -anc {input.painel} \
    -meta {input.meta} \
    -allmod {output.allmod} \
    -figAll {output.figAll} \
    -s {output.summary} \
    -figLat {output.figLat}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule join_syncs:
    input:
        painel_sync = expand(os.path.join(config['analysis_path'], "ancestry/ancestral_chrom_{chrom}.sync"), chrom = config['chrom']),
        sample_sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync")
    output: os.path.join(config['analysis_path'], "ancestry/all_samples.sync")
    params: path_painel_sync = os.path.join(config['analysis_path'], "ancestry")
    shell: "Rscript scripts/R/merge_sync.R -painel {params.path_painel_sync} -sync {input.sample_sync} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule get_all_snp_id:
    input: os.path.join(config['analysis_path'], "ancestry/all_samples.sync")
    output: os.path.join(config['analysis_path'], "ancestry/all_snp_ids.tsv")
    shell: "cut -f 1,2 {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule filter_bs:
    input: 
        ids = os.path.join(config['analysis_path'], "ancestry/all_snp_ids.tsv"),
        bed_sites = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_BS.txt.gz")
    output: os.path.join(config['analysis_path'], "ancestry/all_snp_ids_BS.tsv")
    shell: "zcat {input.bed_sites} | grep -Ff {input.ids} - > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule fix_sample_na:
    input:
        small_bed = os.path.join(config['analysis_path'], "ancestry/all_snp_ids_BS.tsv"),
        all_sync = os.path.join(config['analysis_path'], "ancestry/all_samples.sync")
    output: os.path.join(config['analysis_path'], "ancestry/all_samples_sample_fixed.sync")
    shell: "Rscript scripts/R/fix_sync_sample_na.R -sync {input.all_sync} -bs {input.small_bed} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule filter_file:
    input: os.path.join(config['analysis_path'], "ancestry/all_samples_sample_fixed.sync")
    output: os.path.join(config['analysis_path'], "ancestry/all_snp_ids_filter2.tsv")
    shell: "grep 'NA' {input} | cut -f 1,2 | sed 's/\\t/:/g' | sed '1d' | sed 's/\"//g' > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule make_sync_anc: 
    input:
        ref=config['ref_path'],
        filter_file =  os.path.join(config['analysis_path'], "ancestry/all_snp_ids_filter2.tsv"),
        bams = expand(os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam"), ids = config['haploid_SRAs_dict'].keys())
    output: os.path.join(config['anc_folder'], "call/all_variants_including_pool_samples_sync.sync")
    params: bam_path = os.path.join(config['anc_folder'], "align/"), outdir = os.path.join(config['anc_folder'], "call/")
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf sync \
    --sam-path {params.bam_path} \
    --sam-min-map-qual 20 \
    --sam-min-base-qual 20 \
    --with-header \
    --reference-genome-fasta-file {input.ref} \
    --filter-region-list {input.filter_file} \
    --out-dir {params.outdir} \
    --file-prefix all_variants_including_pool_samples_"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule fix_painel_na:
    input:
        gresync = os.path.join(config['anc_folder'], "call/all_variants_including_pool_samples_sync.sync"),
        allsync = os.path.join(config['analysis_path'], "ancestry/all_samples_sample_fixed.sync"),
        popnames = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/pool_sizes_autosome_with_anc_painel.csv")
    output: allsync = os.path.join(config['analysis_path'], "ancestry/all_samples_sample_and_painel_fixed.sync")
    shell:"Rscript scripts/R/find_non_var_sites_in_painel.R -gre {input.gresync} -all {input.allsync} -pop {input.popnames} -o {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule global_ancetral_pops_FST_autosome:
    input:
        sync = os.path.join(config['analysis_path'], "ancestry/all_samples_sample_and_painel_fixed.sync"),
        pool_sizes = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/pool_sizes_autosome_with_anc_painel.csv"),
        ref = config['ref_path'],
        rename = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/rename_samples_anc_painel.tsv")
    output: os.path.join(config['analysis_path'],"ancestry/FST_ancestry/Genome_FST_autosome_fst-matrix.csv")
    params: outdir=os.path.join(config['analysis_path'], "ancestry/FST_ancestry/"), prefix="Genome_FST_autosome_"
    shell: "/home/vitoria/bin/grenedalf/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta {input.ref} \
    --filter-sample-min-count 1 \
    --filter-sample-min-read-depth 10 \
    --rename-samples-list {input.rename} \
    --filter-region 2L --filter-region 2R --filter-region 3L --filter-region 3R \
    --window-type genome \
    --method unbiased-hudson \
    --window-average-policy valid-loci \
    --pool-sizes {input.pool_sizes} \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}" # \

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule SNP_fst_EU_AFR: #para pegar os snps mais diferenciados entre europa e africa 
    input: 
        sync = os.path.join(config['analysis_path'], "ancestry/ancestral_EU_AFR.sync"),
        ref = config['ref_path'],
        rename = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/rename_samples_EU_AFR.tsv"),
        #regions = "../resources/regions.bed", #*#regions with recombination above 0.5 cM/Mb - from pool et al 2017 - lifted over for ref 6 genome with https://flybase.org/convert/coordinates
        pool_sizes = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/pool_sizes_EU_AFR.csv"),
    output: os.path.join(config['analysis_path'], "ancestry/FST_ancestry/EU_AFR_SNPs_fst.csv")
    params: out_dir = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/")
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --rename-samples-file {input.rename} \
    --window-type single \
    --omit-na-windows \
    --method unbiased-nei \
    --pool-sizes {input.pool_sizes} \
    --verbose \
    --out-dir {params.out_dir} \
    --file-prefix EU_AFR_SNPs_"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule max_fst_list:
    input: 
        painel = os.path.join(config['analysis_path'], "ancestry/raw_ancestry_painel.tsv"),
        fst = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/EU_AFR_SNPs_fst.csv")
    output: os.path.join(config['analysis_path'], "ancestry/FST_ancestry/max_snp_fst_positions.tsv")
    shell: "Rscript scripts/R/max_fst_snp_EU_AFR.R -painel {input.painel} -fst {input.fst} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule highly_diff_SNP_fst: #para pegar os snps mais diferenciados entre europa e africa 
    input: 
        sync = os.path.join(config['analysis_path'], "ancestry/all_samples.sync"),
        pool_sizes = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/pool_sizes_autosome_with_anc_painel.csv"),
        ref = config['ref_path'],
        rename = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/rename_samples_anc_painel.tsv"),
        region = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/max_snp_fst_positions.tsv")
    output: os.path.join(config['analysis_path'], "ancestry/FST_ancestry/highly_diff_SNP_fst.csv")
    params: out_dir = os.path.join(config['analysis_path'], "ancestry/FST_ancestry/")
    shell: "/home/vitoria/bin/grenedalf/bin/grenedalf fst \
    --sync-path {input.sync} \
    --reference-genome-fasta-file {input.ref} \
    --rename-samples-file {input.rename} \
    --filter-region-list {input.region} \
    --window-type single \
    --omit-na-windows \
    --method unbiased-nei \
    --pool-sizes {input.pool_sizes} \
    --verbose \
    --out-dir {params.out_dir} \
    --file-prefix highly_diff_SNP_"
