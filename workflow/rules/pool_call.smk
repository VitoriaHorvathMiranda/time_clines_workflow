IDs = list(grouped_fqs.keys())
IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001")] #deletes samples with low quality from the following analysis

meta_path = config['meta_path']
with open(meta_path, 'r') as file:
    meta = [line.strip().split('\t') for line in file]

list_of_ids = []
for string in IDs_deleted_samples:
    ids = string.split("_")[0]  # Get the part before the first "_"
    list_of_ids.append(ids)

sample_names = []
for string in list_of_ids:
    for row in meta:
        if string == row[0]: #find the row with the id
            sample_names.append(row[1]) #get the pop name
            break  # Stop searching once a match is found
    

#--------------------------------------------------------------------------------------------------------------------------------------------------
# 12 -  - rule mpileup
rule mpileup:
    input:
        fasta=config['ref_path'],
        bam_list_txt = os.path.join(config['align_path'], "bam_list_downsampled.txt")
        #a = os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage/{ids_del}_total_coverage.tsv"),
        #b = "/dados/time_clines/data/seqs/processed/qltctrl/pos_downsample/samtools_depth_stats_include_unmerged_flag.tsv",
        #c = "/dados/time_clines/data/seqs/processed/qltctrl/pre_downsample/samtools_depth_stats_include_unmerged_flag.tsv"
    output: temp(os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup'))
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell:
        "samtools mpileup -A -f {input.fasta} -b {input.bam_list_txt} -q 20 -Q 20 -o {output}"
 #-A flag:Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set
 # precisa colocar essa flag porque os reads das amostras de 2022 não se sobrepõe (eles não estão marcados como "PROPER_PAIR")

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule gzip_mpileup:
    input: os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup')
    output: os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup.gz')
    shell: "gzip {input}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 11 - Variant calling
rule call:
    input: 
        mpileup = os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup.gz'),
        ref = config['ref_path']
    output: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf.gz")
    params:
        outdir= os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15"),
        poolSNP = config['PoolSNP.sh'],
        samples =  ",".join(sample_names)
    threads: 40
    shell: "bash /home/vitoria/bin/PoolSNP/PoolSNP.sh\
    mpileup={input.mpileup}\
    reference={input.ref}\
    names={params.samples} \
    max-cov=0.98\
    min-cov=15\
    min-count=5\
    min-freq=0.001\
    miss-frac=0.1\
    jobs={threads}\
    BS=1\
    base-quality=20\
    allsites=0\
    output={params.outdir}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 13 - vcf unzip
rule unzip_vcf:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf.gz")
    output: temp(os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf"))
    shell: "gunzip -c {input} > {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 14 - get chrom
rule vcf_filter_chrom:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf")
    output: temp(os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_chrom_filtered.vcf"))
    shell: "grep -e $'^2L\t' -e $'^2R\t' -e $'^3L\t' -e $'^3R\t' -e $'^X\t' -e '^#' {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 14 -  #### MARTIN KAPUN SCRIPT (drosEU pipeline - https://github.com/capoony/DrosEU_pipeline)
#identify sites in proximity of InDels with a minimum count of 20 across all samples pooled and mask sites 5bp up- and downstream of InDel
rule Indel_ident:
    input: os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup.gz')
    output: os.path.join(config['align_path'], "InDel-positions_20.txt.gz")
    shell: "python2.7 scripts/python/DetectIndels.py \
--mpileup {input} \
--minimum-count 20 \
--mask 5 \
| gzip > {output}"
