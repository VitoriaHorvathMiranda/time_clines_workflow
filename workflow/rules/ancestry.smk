
#-------------------------------------------------------------------------------------------------------------------------------------------
# the pipeline of the mapping and SNP calling for this vcf is available in the new_samples folder, as is the table detailing the populations used
# 23 central-west African and 22 European (French and Swedish) haploid samples were used 
rule ancestral_painel_prep:
    input: config['ancestral_vcf'],
    output: os.path.join(config['analysis_path'], "ancestry/raw_ancestry_painel.tsv")
    shell: "Rscript scripts/R/raw_ancestral_painel_prep.R -vcf {input} -out {output}"

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
