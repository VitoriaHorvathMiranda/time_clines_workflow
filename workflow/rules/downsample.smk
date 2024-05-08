IDs = list(grouped_fqs.keys())
IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001")] #deletes samples with low quality from the following analysis
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
        dbams = expand(os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam'), ids_del = IDs_deleted_samples),
        dbais = expand(os.path.join(config['align_path'], '{ids_del}_downsample.md.srt.bam.bai'), ids_del = IDs_deleted_samples)
    output: os.path.join(config['align_path'], "bam_list_downsampled.txt")
    wildcard_constraints: ids_del = "|".join(IDs_deleted_samples)
    shell:
        "echo {input.dbams} | sed 's/ /\\n/g' > {output}"

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
    output: "../results/depths_coverage_pre_pos_downsample.jpeg"
    params:
        precp = os.path.join(config['qltctrl_path'], "pre_downsample/samtools_coverage/"),
        poscp = os.path.join(config['qltctrl_path'], "pos_downsample/samtools_coverage/"),
        dpos = config['qltctrl_path']
    shell: "Rscript scripts/R/script_downsample_plot.R -m {input.meta} -precp {params.precp} -poscp {params.poscp} -dpos {params.dpos} -o {output}"