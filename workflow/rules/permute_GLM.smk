PERM_N = [str(i) for i in range(1, 11)]
rule all:
    input:
        expand(os.path.join(config['analysis_path'], "time_GLM_lat/Perm/P_value_perm_n_{perm_n}_{year}.tsv"), perm_n = PERM_N, year = ["97", "0910"]),
        expand(os.path.join(config['analysis_path'], "time_GLM_lat/Perm/enrich_perm_n_{perm_n}_{year}.tsv"), perm_n = PERM_N, year = ["97", "0910"]),
        expand(os.path.join(config['analysis_path'], "time_GLM_lat/Perm/enrich_perm_results_{year}.tsv"), year = ["97", "0910"])

rule make_perm_pvalue:
    input: 
        metadata = config['meta_path'], 
        freqs = os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv"),
        effects = os.path.join(config[''], "time_GLM_lat/effects_q-values_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv")
    output: os.path.join(config['analysis_path'], "time_GLM_lat/Perm/chance_perm_n_{perm_n}_{year}.tsv")
    #group: "permutation"
    wildcard_constraints: perm_n="|".join(PERM_N), year = "(97|0910)"
    shell: "Rscript scripts/R/permuted_lats_glm.R -meta {input.metadata} -freqs {input.freqs} -cyear {wildcards.year} -permN {wildcards.pern_n} -e {input.effects} -out {output}"

rule enrich_perm:
    input:
        p_values = os.path.join(config['analysis_path'], "time_GLM_lat/Perm/P_value_perm_n_{perm_n}_{year}.tsv"),
        vcf_ann = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.ann.vcf")
    output: temp(os.path.join(config['analysis_path'], "time_GLM_lat/Perm/enrich_perm_n_{perm_n}_{year}.tsv"))
    group: "permutation"
    wildcard_constraints: perm_n="|".join(PERM_N), year = "(97|0910)"
    shell: 'echo \"enrichment done\" >> {output}' #"Rscript scripts/R/dummy.R -pvalue {input.p_values} -vcf {input.vcf_ann} -o {output}"

rule analyse_perm:
    input: 
        enrichs = [config['analysis_path'] + "time_GLM_lat/Perm/enrich_perm_n_" + n + "_{year}.tsv" for n in PERM_N]
    output: os.path.join(config['analysis_path'], "time_GLM_lat/Perm/enrich_perm_results_{year}.tsv")
    wildcard_constraints: year = "(97|0910)"
    shell: "cat {input.enrichs} > {output}" #"Rscript scripts/R/dummy2.R -enrichs {input.enrichs} -o {output}"