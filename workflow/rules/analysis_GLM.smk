

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 16 - freq extraction
rule snp_freqs:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: freqs = temp(os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")), snps = os.path.join(config['call_path'], "called_snps_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    shell: "Rscript scripts/R/freq_extraction_pop_ind.R -vcf {input} -snps {output.snps} -o {output.freqs}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 18 - gets metadata and computes NE
rule n_chrom:
    input: 
        metadata = config['meta_path'], 
        depth_file = os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    output: temp(os.path.join(config['call_path'], "NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv"))
    shell: "Rscript scripts/R/n_chrom.R -meta {input.metadata} -df {input.depth_file} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 20 - separates pops based on time
#
rule separate_time_pops:
   input: os.path.join(config['call_path'], "NE_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
   output: expand(os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv"), year=config['YEAR'])
   params:
        path = os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_")

   shell: "Rscript scripts/R/separate_time_pops.R -NE {input} -oPh {params.path}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# separate chom 
# separates files per chrom becase glm script uses too much memory
rule separete_chrom:
    input: os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}.tsv")
    output: os.path.join(config['call_path'], "TimeSEP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_{year}_{chrom}.tsv")
    wildcard_constraints: chrom = "([2-3][LR])|X"
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

