
rule filter_sync:
    input:
        sync = os.path.join(config['anc_folder'], "call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean.sync"),
        ref = config['ref_path'],
        #rename = "../resources/rename_samples.tsv"
    output:
        os.path.join(config['anc_folder'], "call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean_minCount2_minDepth10_totalMAF0.05sync.sync"),
    params: outdir=os.path.join(config['anc_folder'], "call/database_all_zi_eg"), prefix="zi_eg_eu_westafr_all_chrom_biallelic_clean_minCount2_minDepth10_totalMAF0.05"
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "/home/vitoria/bin/grenedalfv6.0/grenedalf-0.6.0/bin/grenedalf sync \
    --sync-path {input.sync} \
    --reference-genome-fasta {input.ref} \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 10 \
    --filter-total-only-snps \
    --filter-total-snp-min-frequency 0.05 \
    --out-dir {params.outdir} \
    --file-prefix {params.prefix}" 
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule anc_pca:
    input: 
        vcf = os.path.join(config['anc_folder'], "call/database_all_zi_eg/zi_eg_eu_westafr_{chrom}_biallelic_clean.vcf"),
        sync = os.path.join(config['anc_folder'], "call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean.sync"),
        inv = "../resources/inversions_anc_samples.csv"
    output: "../results/anc_pca_{chrom}.pdf"
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "Rscript scripts/R/anc_pca.R -vcf {input.vcf} -sync {input.sync} -chrom {wildcards.chrom} -inv {input.inv} -pca {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule join_sync_f4:
    input: 
        anc_sync = os.path.join(config['anc_folder'], "call/database_all_zi_eg/zi_eg_eu_westafr_all_chrom_biallelic_clean.sync"),
        my_sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: os.path.join(config['analysis_path'], "ancestry/f4ratio/{chrom}.sync")
    params: path = os.path.join(config['analysis_path'], "ancestry/f4ratio/")
    shell: "Rscript scripts/R/join_all_sync.R -vcf {input.vcf} -ancSync {input.anc_sync} -mySync {input.my_sync} -o {output} -opath {params.path}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule make_anc_freq_table:
    input:
        sync = os.path.join(config['analysis_path'], "ancestry/f4ratio/all.sync"),
        ref = config['ref_path'],
        rename = os.path.join(config['analysis_path'], "ancestry/f4ratio/rename_all.tsv")
    output:
        os.path.join(config['analysis_path'], "ancestry/f4ratio/all_freq_table_frequency.csv"),
    params: outdir=os.path.join(config['analysis_path'], "ancestry/f4ratio"), prefix="all_freq_table_"
    #wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "/home/vitoria/bin/grenedalfv6.0/grenedalf-0.6.0/bin/grenedalf frequency \
    --sync-path {input.sync} \
    --reference-genome-fasta {input.ref} \
    --rename-samples-list {input.rename} \
    --out-dir {params.outdir} \
    --write-sample-alt-freq \
    --write-sample-read-depth \
    --file-prefix {params.prefix}" 

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule f4ratio:
    input: 
        sync = os.path.join(config['analysis_path'], "ancestry/f4ratio/{chrom}.sync"),
        sizes = "../resources/pool_sizes_autosome.csv",
        meta = config['meta_path'],
    output: 
        o = "../results/f4ratio_estimates_{chrom}.tsv", 
        plot = "../results/f4ratio_per_latitude_{chrom}.pdf"
    wildcard_constraints: chrom = config['chrom'].append("all")
    shell: "Rscript scripts/R/f4ratio.R -sync {input.sync} -sizes {input.sizes} -meta {input.meta} -o {output.o} -plot {output.plot}"

