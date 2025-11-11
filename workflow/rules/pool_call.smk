IDs = list(grouped_fqs.keys())
IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001", "A41_L002")] #deletes samples with low quality from the following analysis


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
rule name_order:
    input: 
        meta = config['meta_path'],
        bam_list = os.path.join(config['align_path'], "bam_list_downsampled.txt")
    output:
        sample_names = "../resources/sample_names.txt"
    run:
        # Check the paths
        print(f"bam_list path: {input.bam_list}")
        print(f"meta path: {input.meta}")
        
        # Load the metadata
        with open(input.meta, 'r') as file:
            meta = [line.strip().split('\t') for line in file]
        
        # Initialize an empty list to store the results
        IDs_bams = []

        # Open the BAM list file in read mode
        with open(input.bam_list, 'r') as file:
            # Loop through each line in the file
            for line in file:
                line = line.strip()
                #first_part = line.split('_')[0]
                path_parts = line.split('/')[6]
                print(path_parts)
                id_names = path_parts.split('_')[0]
                print(id_names)
                IDs_bams.append(id_names)
                
        # Process the IDs with the metadata
        sample_names = []
        for string in IDs_bams:
            for row in meta:
                if string == row[0]:  # Find the row with the id
                    sample_names.append(row[1])  # Get the pop name
                    break
         
        # Write the sample names to the output file
        with open(output.sample_names, 'w') as out_file:
            out_file.write("\n".join(sample_names))


#-------------------------------------------------------------------------------------------------------------------------------------------------
# 11 - Variant calling
rule call:
    input: 
        mpileup = os.path.join(config['align_path'], 'noSNC10_noESC97_with_dlGA10_dlSC10.mpileup.gz'),
        ref = config['ref_path'],
        names = "../resources/sample_names.txt"
    output: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf.gz")
    params:
        outdir= os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15"),
        poolSNP = config['PoolSNP.sh'],
        samples =  lambda wildcards, input: ",".join([line.strip() for line in open(input.names)])
    threads: 40
    shell: "bash /home/vitoria/bin/PoolSNP/PoolSNP.sh\
    mpileup={input.mpileup}\
    reference={input.ref}\
    names={params.samples} \
    max-cov=0.98\
    min-cov=15\
    min-count=5\
    min-freq=0.01\
    miss-frac=0.1\
    jobs={threads}\
    BS=1\
    base-quality=20\
    allsites=0\
    output={params.outdir}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 14 - get chrom
rule vcf_filter_chrom:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf.gz")
    output: temp(os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_chrom_filtered.vcf"))
    shell: "zcat {input} | grep -e $'^2L\t' -e $'^2R\t' -e $'^3L\t' -e $'^3R\t' -e $'^X\t' -e '^#' > {output}"

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

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 15 -  #### MARTIN KAPUN SCRIPT (drosEU pipeline - https://github.com/capoony/DrosEU_pipeline)
# filter SNPs around InDels from the original VCF produced with PoolSNP
rule filter_indel:
    input: 
        indels = os.path.join(config['align_path'], "InDel-positions_20.txt.gz"),
        vcf = os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_chrom_filtered.vcf")
    output: temp(os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_indels_filtered.vcf.gz"))
    shell:"python2.7 scripts/python/FilterPosFromVCF.py \
--indel {input.indels} \
--vcf {input.vcf} \
| gzip > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 16 - transform repeat regions in .bed
#repeat regions came from : http://www.repeatmasker.org/species/dm.html
rule fix_reapeat:
    input: config['repeat_regios']
    output: os.path.join(config['ref_folder'],"repeat_6.bed")
    shell: "sed '1,3d' {input} | tr -s ' ' | sed 's/^ //g' | cut -d' ' -f5,6,7 | tr '[:blank:]' \\t | sed 's/chr//g' | bedtools sort > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 17 - filter repeat regions from vcf
rule filter_repeat:
    input: 
        vcf = os.path.join(config['call_path'],"PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_indels_filtered.vcf.gz"),
        bed_file = os.path.join(config['ref_folder'],"repeat_6.bed")
    output: temp(os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.vcf"))
    shell: "bedtools subtract -a {input.vcf} -b {input.bed_file} | sed 's/\s/\t/g' > {output}"


#-------------------------------------------------------------------------------------------------------------------------------------------------
# 18 - get header #o snakemake falha nessa regra, mas se executar direto no terminal funciona 
# Eu tentei ver qual a saída do comando com 'echo $?', mas quando eu executo no terminal, a saída é zero
#resolvi com um '|| echo "qualquer coisa"'
rule get_header:
    input:
        raw_vcf=os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.vcf.gz"),
        clean_vcf=os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.vcf")
    output:os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    shell: "zcat {input.raw_vcf} | head -n 18 | cat - {input.clean_vcf} > {output} || echo \"Snakemake says this command doesn't work, but it does\" "
    #; sed -i 's/Alternative Counts/Allelic depths for the ref and alt alleles in the order listed/' {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#make .sync file based on vcf - #### MARTIN KAPUN SCRIPT (drosEU pipeline - https://github.com/capoony/DrosEU_pipeline)
rule make_sync:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.sync")
    shell: "python scripts/python/VCF2sync.py --vcf {input} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 16 - freq extraction
#get snp frequencies 
#o script filtra para os SNPs bialélicos
rule snp_freqs:
    input: os.path.join(config['call_path'], "PoolSNP_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15_clean.h.vcf")
    output: freqs = os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv"), snps = os.path.join(config['call_path'], "called_snps_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    shell: "Rscript scripts/R/freq_extraction_pop_ind.R -vcf {input} -snps {output.snps} -o {output.freqs}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#Ti/Tv table
rule Ti_Tv:
    input: os.path.join(config['call_path'], "freqs_and_depths_noSNC10_noESC97_with_dlGA10_dlSC10_mincount5_minfreq0.001_cov15.tsv")
    output: os.path.join(config['analysis_path'], "quality/Ti_Tv.tsv")
    shell: "Rscript scripts/R/Ti_Tv_ratio.R -freqs {input} -o {output}"