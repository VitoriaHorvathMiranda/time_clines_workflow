IDs = list(grouped_fqs.keys())
IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001", "A41_L002")] #deletes samples with low quality from the following analysis
male_labels = [id + "_L001" for id in config['male_ids']]
IDs_female_pools = [e for e in IDs_deleted_samples if e not in male_labels]
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 11 - rule bam.txt
rule bam_txt:
    input: 
        bams = expand(os.path.join(config['align_path'], '{ids_del}.md.srt.bam'), ids_del = IDs_deleted_samples), 
        bais = expand(os.path.join(config['align_path'], '{ids_del}.md.srt.bam.bai'),ids_del = IDs_deleted_samples),
    output: os.path.join(config['align_path'], "bam_list.txt")
    #wildcard_constraints: id="[a-zA-Z0-9]+", lane="\d"
    shell:
        "echo {input.bams} | sed 's/ /\\n/g' > {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
# - computes depth per position
rule compute_depth_stats:
    input: bam_list_txt= os.path.join(config['align_path'], "bam_list.txt"), 
    output: os.path.join(config['qltctrl_path'], "pre_downsample/samtools_depth_stats_include_unmerged_flag.tsv")
    #wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell:
        "samtools depth -a -H -f {input.bam_list_txt} -q 20 -Q 20 -o {output}"
#-------------------------------------------------------------------------------------------------------------------------------------------------
# compute total coverage per chromosome, mean depth and mean quality
rule total_coverage:
    input: 
        bams = os.path.join(config['align_path'], '{ids_del}.md.srt.bam'), 
        bais = os.path.join(config['align_path'], '{ids_del}.md.srt.bam.bai')
    output: os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage/{ids_del}_total_coverage.tsv")
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell:"samtools coverage -q 20 -Q 20 -d 200 -o {output} {input.bams}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#samtools idxstats – reports alignment summary statistics
#get the number of mapped and unmapped reads
rule align_stats:
    input: 
        bams = os.path.join(config['align_path'], '{ids_del}.md.srt.bam'), 
        bais = os.path.join(config['align_path'], '{ids_del}.md.srt.bam.bai')
    output: os.path.join(config['qltctrl_path'], "pre_downsample/samtools_idxstats/{ids_del}_alignStats.tsv")
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell: "samtools idxstats {input.bams} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#writes a tsv file with the downsample fraction need to downsample reads to a mean depth of 25 
rule downsample_fraq:
    input:
        meta = config['meta_path'],
        stats_path = os.path.join(config['qltctrl_path'], 'pre_downsample/samtools_idxstats/'),
        cov_path = os.path.join(config['qltctrl_path'], 'pre_downsample/samtools_coverage/'),
        stats = expand(os.path.join(config['qltctrl_path'], "pre_downsample/samtools_idxstats/{ids_del}_alignStats.tsv"), ids_del = IDs_deleted_samples), 
        cov = expand(os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage/{ids_del}_total_coverage.tsv"), ids_del = IDs_deleted_samples) 
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    output: os.path.join(config['qltctrl_path'], "downsample_info.tsv")
    shell: "Rscript scripts/R/script_reads_cutoff.R -meta {input.meta} -stats {input.stats_path} -cov {input.cov_path} -o {output}"#
    # the script computs the read fraq we need so that the mean depth of the reads is equal to 25X
    # the script takes all files with specific string in the input folders, so if there are samples that aren't part of the ids_del list, it will also compute the read fraq for them

#--------------------------------------------------------------------------------------------------------------------------------------------------
rule downsample_samtools:
    input: bam = os.path.join(config['align_path'], '{ids_del}.md.srt.bam'), frac_stats = os.path.join(config['qltctrl_path'], "downsample_info.tsv")
    output: os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam')
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    params: population = lambda wildcards: wildcards.ids_del.split("_")[0]
    threads: 5
    shell: "frac_reads=$(cat {input.frac_stats} | awk '$1 == \"{params.population}\" {{print ($9 >= 1) ? 1 : $9}}') && \
    samtools view -@{threads} -bs $frac_reads {input.bam} > {output} || cp {input.bam} {output}"
    #($9 >= 1) is a condition that checks if the value of the 9th field is greater than or equal to 1
    #? 1 : $9 -> if the condition (9th field >= 1) is true, it evaluates to 1. If the condition is false, it evaluates to the value of the 9th field itself.
    # || cp -> the command will fail when $frac_reads == 1, then it will only copy the original bam file (which has fewer reads than the threshold)

#--------------------------------------------------------------------------------------------------------------------------------------------------
rule index_downsampled:
    input: os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam')
    output: os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam.bai')
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    threads: 4
    shell:
        "samtools index -@{threads} {input}"

#--------------------------------------------------------------------------------------------------------------------------------------------------
# 11 - rule bam.txt
rule bam_txt_downsampled:
    input: 
        dbams = expand(os.path.join(config['align_path'], '{ids_fe}_downsample.md.srt.bam'), ids_fe = IDs_female_pools),
        dbams_males = expand(os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam'), ids_male = config['male_ids']),
        dbais = expand(os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam.bai'), ids_del = IDs_deleted_samples),
        dbais_males = expand(os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam.bai'), ids_male = config['male_ids'])
    output: os.path.join(config['align_path'], "bam_list_downsampled.txt")
    wildcard_constraints: ids_fe = "|".join(IDs_female_pools), ids_male = "|".join(config['male_ids'])
    shell:
        "echo {input.dbams} {input.dbams_males} | sed 's/ /\\n/g' > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# - computes depth per position pos downsampling
rule compute_depth_stats_downsampled:
    input: dbam_list_txt= os.path.join(config['align_path'], "bam_list_downsampled.txt")
    output: os.path.join(config['qltctrl_path'], "pos_downsample/samtools_depth_stats_include_unmerged_flag.tsv" ) 
    resources:
        mem_mb=34000
    shell:
        "samtools depth -a -H -f {input.dbam_list_txt} -q 20 -Q 20 -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# compute total coverage
rule total_coverage_downsampled:
    input: 
        dbams = os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam'),
        dbais = os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam.bai')
    output: os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage/{ids_del}_total_coverage.tsv")
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell:"samtools coverage -q 20 -Q 20 -d 200 -o {output} {input.dbams}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule get_male_X_bam:
    input: 
        bams = os.path.join(config['align_path'], '{ids_male}_L001.md.srt.bam'), 
        bais = os.path.join(config['align_path'], '{ids_male}_L001.md.srt.bam.bai')
    output: temp(os.path.join(config['align_path'], '{ids_male}_only_chrom_X.md.srt.bam')), 
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    shell: "samtools view -b {input.bams} X > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule down_fraq_male_X_bam:
    input: cov = expand(os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage/{male_ids}_L001_total_coverage.tsv"), male_ids = config['male_ids'])
    output:  os.path.join(config['qltctrl_path'], "downsample_info_chromX.tsv")
    params: 
        cov_path = os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage"),
        pops= "_".join(config['male_ids'])
    shell: "Rscript scripts/R/downsample_fraq_male_pools.R -cov {params.cov_path} -mpools {params.pops} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule downsample_X_chrom:
    input: 
        xbam = os.path.join(config['align_path'], '{ids_male}_only_chrom_X.md.srt.bam'),
        frac = os.path.join(config['qltctrl_path'], "downsample_info_chromX.tsv")
    output: temp(os.path.join(config['align_path'], '{ids_male}_only_chrom_X_downsampled.md.srt.bam'))
    params: pops = "|".join(config['male_ids'])
    threads: 4
    shell: "frac_reads=$(cat {input.frac} | awk '$10 == \"{wildcards.ids_male}\" {{print ($11 >= 1) ? 1 : $11}}') && \
    samtools view -@{threads} -bs $frac_reads {input.xbam} > {output} || cp {input.xbam} {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule delete_X_male_down:
    input: 
        dbams = os.path.join(config['align_path'], '{ids_male}_L001_downsample.md.srt.bam'),
        X_filter = "../resources/CHORM_X.bed",
    output: 
        noX = temp(os.path.join(config['align_path'], '{ids_male}_noX_downsample.md.srt.bam'))
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    threads: 4
    shell: "samtools view -@{threads} -b {input.dbams} -L {input.X_filter} -U {output.noX}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule merge_male_bams:
    input: 
        xbam = os.path.join(config['align_path'], '{ids_male}_only_chrom_X_downsampled.md.srt.bam'),
        noXbam = os.path.join(config['align_path'], '{ids_male}_noX_downsample.md.srt.bam')
    output: temp(os.path.join(config['align_path'], '{ids_male}_all_male_chrom.md.srt.bam'))
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    threads: 4
    shell: "samtools cat {input.noXbam} {input.xbam} -o {output}"
    #"java -jar /home/vitoria/bin/picard.jar MergeSamFiles I={input.noXbam} I={input.xbam} O={output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule male_sort:
    input: os.path.join(config['align_path'], '{ids_male}_all_male_chrom.md.srt.bam')
    output: os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam')
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    shell: "samtools sort -o {output} {input}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule male_index:
    input: os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam')
    output: os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam.bai')
    threads: 5
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    shell: "samtools index -@{threads} {input}"


#-------------------------------------------------------------------------------------------------------------------------------------------------
rule male_covstats:
    input: os.path.join(config['align_path'], '{ids_male}_all_male_chrom_sorted.md.srt.bam')
    output: os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage/{ids_male}_all_chrom_total_coverage.tsv")
    wildcard_constraints: ids_male = "|".join(config['male_ids'])
    shell: "samtools coverage -q 20 -Q 20 -d 200 -o {output} {input}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#compute depth per 80kb windows
rule depth_per_window:
    input:
        pre = os.path.join(config['qltctrl_path'], "pre_downsample/samtools_depth_stats_include_unmerged_flag.tsv" ),
        pos = os.path.join(config['qltctrl_path'], "pos_downsample/samtools_depth_stats_include_unmerged_flag.tsv" )
    output: os.path.join(config['qltctrl_path'], "depth_per_80kb_window_{chrom}.tsv")
    wildcard_constraints: chrom="|".join(config['chrom'])
    shell: "Rscript scripts/R/script_depth_per_pos.R -predf {input.pre} -posdf {input.pos} -c {wildcards.chrom} -o {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
#plot depth per position before and after downsample and mean depth per coverage
rule plot_depth:
    input:
        meta = config['meta_path'],
        pre = expand(os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage/{ids_del}_total_coverage.tsv"), ids_del = IDs_deleted_samples),
        pos = expand(os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage/{ids_del}_total_coverage.tsv"), ids_del = IDs_deleted_samples),
        depths = expand(os.path.join(config['qltctrl_path'], "depth_per_80kb_window_{chrom}.tsv"), chrom = config['chrom'])
    output: 
        big_plot = "../results/depths_coverage_pre_pos_downsample.jpeg",
        boxplots = "../results/mean_depths_boxplot.jpeg"
    params:
        precp = os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage"),
        poscp = os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage"),
        dpos = config['qltctrl_path'],
        pops= "_".join(config['male_ids'])
    shell: "Rscript scripts/R/script_downsample_plot.R -m {input.meta} -precp {params.precp} -poscp {params.poscp} -mIDs {params.pops} -dpos {params.dpos} -o {output.big_plot} -box {output.boxplots}"