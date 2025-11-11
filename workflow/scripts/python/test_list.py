import os
import yaml

with open("../../../config/config.yaml") as f:
    config = yaml.safe_load(f)

w_fst_cutoff_merged = []
for chrom in config["chrom"]:
    for local in config["local_to_pairs"]:
        for pair in config["local_to_pairs"][local]:
            filename = 'fst/window/merged_cutoffs_0.001_' + local + '_100_' + pair + '_' + chrom + '.tsv'
            full_path = os.path.join(config['analysis_path'], filename)
            w_fst_cutoff_merged.append(full_path)

#for path in w_fst_cutoff_merged:
#    print(path)

#_all_chroms_with_genes.tsv
w_fst_cutoff_merged_all_chrom = []
for local in config["local_to_pairs"]:
    for pair in config["local_to_pairs"][local]:
        filename = 'fst/window/merged_cutoffs_0.001_' + local + '_100_' + pair + '_all_chroms_with_genes.tsv'
        full_path = os.path.join(config['analysis_path'], filename)
        w_fst_cutoff_merged_all_chrom.append(full_path)

window_size = []
for local in config["local_to_pairs"]:
    for chrom in config["chrom"]:
        filename = 'fst/window/BP_windown_size_' + local + '_100_' + chrom + '.tsv'
        full_path = os.path.join(config['analysis_path'], filename)
        window_size.append(full_path)

for path in window_size:
    print(path)
