
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
rule get_inv_markers:
    input:
        pos = "../resources/inversion_markers_v6.txt_pos", #file from https://github.com/capoony/DrosEU_pipeline
        sync = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync")
    output: os.path.join(config['analysis_path'], "inversions/inversion_markers.sync")
    shell: "python2 scripts/python/OverlapSNPs.py --source {input.pos} --target {input.sync} > {output}" #MARTIN KAPUN SCRIPT from https://github.com/capoony/DrosEU_pipeline

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#makes a file with pop name and corresponding pop column in .sync file
rule sync_col:
    input: "../resources/rename_samples.tsv"
    output: "../resources/sync_pop_position.txt" #not a txt
    shell: "awk '{{print $2 \"\\t\" NR+2}}' {input} | sed '1i\\\\' > {output}" 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
rule Inv_freq:
    input: 
        inv_sync = os.path.join(config['analysis_path'], "inversions/inversion_markers.sync"),
        sync_names = "../resources/sync_pop_position.txt",
        inv = "../resources/inversion_markers_v6.txt" #file from https://github.com/capoony/DrosEU_pipeline
    output: os.path.join(config['analysis_path'], "inversions/inversion_freq.tsv")
    shell: "python scripts/python/InvFreq.py --sync {input.inv_sync} --meta {input.sync_names} --inv {input.inv} > {output}"

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
rule plot_invFreq:
    input: 
        meta = config['meta_path'],
        inv_freq = os.path.join(config['analysis_path'], "inversions/inversion_freq.tsv")
    output: plot = "../results/Inv_freq_per_lat.jpeg", glms = "../results/Inv_freq_per_lat.txt"
    shell: "Rscript scripts/R/Inv_freq_plot.R -meta {input.meta} -invf {input.inv_freq} -out {output.plot} -glm {output.glms}"