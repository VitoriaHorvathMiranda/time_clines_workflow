

PERM_N = [str(i) for i in range(1, 5041)] #a primeira permutação são os latutudes certas
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 18 - gets metadata and computes NE
rule n_chrom:
    input: 
        metadata = config['meta_path'], 
        depth_file = os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    output: os.path.join(config['call_path'], "NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    shell: "Rscript scripts/R/n_chrom.R -meta {input.metadata} -df {input.depth_file} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 20 - separates pops based on time
#
rule separate_time_pops:
   input: os.path.join(config['call_path'], "NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
   output: expand(os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv"), year=config['CLINAL_YEAR'])
   wildcard_constraints: year= "|".join(config['CLINAL_YEAR'])
   params:
        path =  os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_")

   shell: "Rscript scripts/R/separate_time_pops.R -NE {input} -oPh {params.path}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# separate chom 
# separates files per chrom becase glm script uses too much memory
rule separete_chrom:
    input: os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv")
    output: temp(os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv"))
    wildcard_constraints: chrom = "([2-3][LR])|X", year= "|".join(config['CLINAL_YEAR'])
    shell:"awk '($4==\"{wildcards.chrom}\")' {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 21 - GLM
#gets p-values
rule glm:
    input: os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv")
    output: os.path.join(config['analysis_path'], "time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv")
    wildcard_constraints: chrom = "([2-3][LR])|X", year= "|".join(config['CLINAL_YEAR'])
    resources:
        mem_mb= 20000
    shell: "Rscript scripts/R/glm_script.R --timePops {input} --output {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#
rule make_glm_plots:
    input: expand(os.path.join(config['analysis_path'], "time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv"), year=config['CLINAL_YEAR'], chrom = config['chrom'])
    output: 
        m97 = "../results/manhattan_97.jpeg",
        m0910 = "../results/manhattan_0910.jpeg",
        mfl2 = "../results/manhattan_0910_FL2.jpeg",
        ht = "../results/P-values_hist.jpeg",
        q97 = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv"),
        q0910 = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910.tsv"),
        q0910FL2 = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910FL2.tsv"),
        qhists = "../results/hist_p-values_withPi.jpeg",
        qPlot97 = "../results/q_Plot97.jpeg",
        qPlot0910 = "../results/q_Plot0910.jpeg",
        qPlot0910FL2 = "../results/q_Plot0910FL2.jpeg"
    params: 
        path= os.path.join(config['analysis_path'], "time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_"),
        Opath= os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_")
    wildcard_constraints: chrom = "([2-3][LR])|X", year= "|".join(config['CLINAL_YEAR']), 
    shell: "Rscript scripts/R/p_value_hist_manhattan.R \
    -path {params.path} \
    -Opath {params.Opath} \
    -m97 {output.m97} -m0910 {output.m0910} -m0910FL2 {output.mfl2} \
    -hist {output.ht}\
    -qhist {output.qhists}\
    -qp97 {output.qPlot97} -qp0910 {output.qPlot0910} -qp0910FL2 {output.qPlot0910FL2}"

#-------------------------------------------------------------------------------------------------------------------------------------------------

#rule get_qvalues:
#    input: expand(os.path.join(config['analysis_path'], "time_GLM_lat/p-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv"), year=config['CLINAL_YEAR'], chrom = config['chrom'])
#    output:
#        q97 = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv"),
#        q0910 = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910.tsv"),
#        hists = "../results/hist_p-values_withPi.jpeg",
#        qPlot97 = "../results/q_Plot97.jpeg",
#        qPlot0910 = "../results/q_Plot0910.jpeg"
#    params: path_results = os.path.join(config['analysis_path'], "time_GLM_lat")
#    shell: "Rscript scripts/R/q_values_script.R -snps {params.path_results} -o97 {output.q97} -o0910 {output.q0910} -hp {output.hists} -qp97 {output.qPlot97} -qp0910 {output.qPlot0910}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule annotate:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.ann.vcf")
    shell: "java -Xmx8g -jar /home/vitoria/bin/snpEff/snpEff.jar -ud 1000 BDGP6.32.105 {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule get_effects: #get the strongest effect of each snp (the first from SNPeff)
    input:
        annvcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.ann.vcf"),
        effects = "../resources/effects.tsv",
        q = os.path.join(config['analysis_path'], "time_GLM_lat/q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
    output: os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR'])
    shell: "Rscript scripts/R/get_effects.R -vcf {input.annvcf} -ef {input.effects} -qv {input.q} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule genomic_region_odr:
    input: 
        qvalues = os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv"),
    output: os.path.join(config['analysis_path'], "time_GLM_lat/genomic_region_odr/odds_ratio_{clinal_year}_clinal_at_{fdr}.tsv")
    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR']), fdr = "|".join(config['FDR_cutoffs'])
    shell: "Rscript scripts/R/odds_ratio_to_random_sampling.R  -qv {input.qvalues} -cut {wildcards.fdr} -out {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule plot_genomic_region_odr:
    input: 
        odr97 = os.path.join(config['analysis_path'], "time_GLM_lat/genomic_region_odr/odds_ratio_97_clinal_at_0.1.tsv"),
        odr0910 = os.path.join(config['analysis_path'], "time_GLM_lat/genomic_region_odr/odds_ratio_0910_clinal_at_0.1.tsv"),
    output: "../results/odr_genomic_region_random_sampling.jpeg"
    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR']), 
    shell: "Rscript scripts/R/plot_odd_ratio_sampling.R  -odr97 {input.odr97} -odr0910 {input.odr0910} -out {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule join_chances:
#    input: [config['analysis_path'] + "time_GLM_lat/Perm/chance_perm_n_" + n + "_{year}.tsv" for n in PERM_N]
#    output: os.path.join(config['analysis_path'], "time_GLM_lat/Perm/joined_perm_chances_{year}.tsv")
#    params: path = os.path.join(config['analysis_path'], "time_GLM_lat/Perm")
#    shell: "Rscript scripts/R/join_perm_chances.R -cyear {wildcards.year} -pc {params.path}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule enrichment_odds_ratio:
    input: os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
    output:
        ench= "../results/odds_ratio_{clinal_year}.jpeg",
        table_odds = os.path.join(config['analysis_path'], "time_GLM_lat/odds_ratio_{clinal_year}.tsv"),
        bar_plot = "../results/effect_summary_{clinal_year}.jpeg"
    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR'])
    shell: "Rscript scripts/R/enrichment_effects.R -eff {input} -pen {output.ench} -out {output.table_odds} -cp {output.bar_plot}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
rule enrichment_to_integenic:
    input: os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
    output:
        ench= "../results/odds_ratio_to_integenic_{clinal_year}.jpeg",
        table_odds = os.path.join(config['analysis_path'], "time_GLM_lat/odds_ratio_to_intergenic_{clinal_year}.tsv")
    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR'])
    shell: "Rscript scripts/R/enrichment_to_intergenic.R -eff {input} -pen {output.ench} -out {output.table_odds}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule make_perm_pvalue:
#    input: 
#        metadata = config['meta_path'], 
#        freqs = os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv"),
#        effects = os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv")
#    output: os.path.join(config['analysis_path'], "time_GLM_lat/Perm/chance_perm_n_{perm_n}_{year}.tsv")
#    #group: "permutation"
#    wildcard_constraints: perm_n="|".join(PERM_N), year = "(97|0910)"
#    shell: "Rscript scripts/R/permuted_lats_glm.R -meta {input.metadata} -freqs {input.freqs} -cyear {wildcards.year} -permN {wildcards.perm_n} -e {input.effects} -out {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule join_chances:
#    input: [config['analysis_path'] + "time_GLM_lat/Perm/chance_perm_n_" + n + "_{year}.tsv" for n in PERM_N]
#    output: os.path.join(config['analysis_path'], "time_GLM_lat/Perm/joined_perm_chances_{year}.tsv")
#    params: path = os.path.join(config['analysis_path'], "time_GLM_lat/Perm")
#    shell: "Rscript scripts/R/join_perm_chances.R -cyear {wildcards.year} -pc {params.path}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule enrichment_odds_ratio:
#    input: os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
#    output:
#        ench= "../results/odds_ratio_{clinal_year}.jpeg",
#        table_odds = os.path.join(config['analysis_path'], "time_GLM_lat/odds_ratio_{clinal_year}.tsv"),
#        bar_plot = "../results/effect_summary_{clinal_year}.jpeg"
#    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR'])
#    shell: "Rscript scripts/R/enrichment_effects.R -eff {input} -pen {output.ench} -out {output.table_odds} -cp {output.bar_plot}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule enrichment_to_integenic:
#    input: os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{clinal_year}.tsv")
#    output:
#        ench= "../results/odds_ratio_to_integenic_{clinal_year}.jpeg",
#        table_odds = os.path.join(config['analysis_path'], "time_GLM_lat/odds_ratio_to_intergenic_{clinal_year}.tsv")
#    wildcard_constraints: clinal_year = "|".join(config['CLINAL_YEAR'])
#    shell: "Rscript scripts/R/enrichment_to_intergenic.R -eff {input} -pen {output.ench} -out {output.table_odds}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule get_clinal_candidate_snps:
#    input: 
#        e97 = os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_97.tsv"),
#        e0910 = os.path.join(config['analysis_path'], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_0910.tsv"),
#    output: 
#        only97= os.path.join(config['analysis_path'], "time_GLM_lat/GO_analysis/only_clinal_in_97_at_{fdr}_fdr_candidate_snps.tsv"),
#        only0910 = os.path.join(config['analysis_path'], "time_GLM_lat/GO_analysis/only_clinal_in_0910_at_{fdr}_fdr_candidate_snps.tsv")
#    wildcard_constraints: fdr = "|".join(config['FDR_cutoffs'])
#    shell: "Rscript scripts/R/get_clinal_candidateSNPs.R -e97 {input.e97} -e0910 {input.e0910} -fdr {wildcards.fdr} -clinal97 {output.only97} -clinal0910 {output.only0910}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# make go assossiation file
#rule GO_association:
#    input: gtf=os.path.join(config['ref_folder'], "dmel-6.55.gtf")
#    output: os.path.join(config['ref_folder'],"GO_gene_association.tsv")
#    shell: "Rscript scripts/R/GO_association.R -gtf {input.gtf} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#
#rule GO_analysis_clinal_s:
#    input:
#        gtf=os.path.join(config['ref_folder'], "dmel-6.55.gtf"),
#        goassociation =  os.path.join(config['ref_folder'],"GO_gene_association.tsv"),
#        total_snps = os.path.join(config['call_path'],"PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync"),
#        cand_snps= os.path.join(config['analysis_path'], "time_GLM_lat/GO_analysis/only_clinal_in_{year}_at_{fdr}_fdr_candidate_snps.tsv")
#    output: os.path.join(config['analysis_path'], "time_GLM_lat/GO_analysis/GO_enrichment_mode_{mode}_gendef_{gen_def}_only_clinal_{year}_at_{fdr}_fdr.tsv")
#    threads: 8
#    params: output_name="GO_enrichment_mode_{mode}_gendef_{gen_def}_only_clinal_{year}_at_{fdr}_fdr.tsv"
#    wildcard_constraints: year = "|".join(config['CLINAL_YEAR']), fdr = "|".join(config['FDR_cutoffs']), mode = "(gene|snp)", gen_def = "updownstream1000"
#    shell: "java -Xmx4g -jar /home/vitoria/bin/Gowinda-1.12.jar --snp-file {input.total_snps} \
#    --candidate-snp-file {input.cand_snps} \
#    --gene-set-file {input.goassociation} \
#    --annotation-file {input.gtf} \
#    --simulations 100000 \
#    --min-significance 1 \
#    --gene-definition {wildcards.gen_def} \
#    --threads {threads} \
#    --output-file {params.output_name} \
#    --mode {wildcards.mode} \
#    --min-genes 5 || mv {params.output_name} {output}"
