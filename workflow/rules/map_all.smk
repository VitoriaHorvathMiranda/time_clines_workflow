meta = pd.read_csv(config['meta_path'], sep="\t")
read_groups = {}
sample_names = {}
def header_to_readgroup(sample, header):
    rid = ".".join(str(header).split(":")[2:4])
    return f'\"@RG\\tID:{rid}\\tSM:{sample}\\tPL:Illumina\"'
for fqpath, fqname in zip(r1_fqs, fqs_pref): # zip => exemplo: a = (1, 2, 3, 4); d = (12, 13); zip(a, d) = [(1, 12), (2, 13)]
    sampleid = fqname.split('_')[0]
    sid = sampleid
    if sampleid == "09":
        sid="9"
    samplename = meta[meta.seq_label==sid].population.values[0]
    sample_names[sampleid] = samplename
    f = gzip.open(fqpath, 'r')
    line = f.readline()
    read_groups[fqname] = header_to_readgroup(samplename, line)
    f.close()

#----------------------------------------------------------------------------------------------------------------------------------------------
rule ungz:
    input: r=os.path.join(config['processed_path'], "{ids}_merged.fastq.gz")
    output: temp(os.path.join(config['processed_path'], "{ids}_merged.fastq"))
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
    threads: 4
    shell:
        "pigz -f -p {threads} -d {input} > {output}"
        # Compress or expand files
            # -f (--force) = Force overwrite, compress .gz, links, and to terminal.
            # -p (--processes) = Allow up to n processes (default is the number of online processors).
            # -d (--decompress or --uncompress) = Decompress the compressed input.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 3.0 - Regra eh utilizada no input da rule bwa
rule bwaindex:
    input: config['ref_path']
    output: config['ref_index_path']
    shell:
        "bwa-mem2 index {input}; " # Index database sequences in the FASTA format.
        "samtools faidx {input}; " # Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence.
                                 # If no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk.
                                 # If regions are specified, the subsequences will be retrieved and printed to stdout in the FASTA format.
        "gatk CreateSequenceDictionary -R {output}" # Create a SAM/BAM file from a fasta containing reference sequence.
                                                    # The output SAM file contains a header but no SAMRecords, and the header
                                                    # contains only sequence records.

# 3.1 - Mapeia os reads contra a referencia - cria os arquivos .sam
rule bwa:
    input: r=os.path.join(config['processed_path'], "{ids}_merged.fastq"), refi=config['ref_index_path'], ref=config['ref_path']
    params: rg = lambda wildcards: read_groups[wildcards.ids]
    output: temp(os.path.join(config['align_path'], "{ids}.sam"))
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
    threads: 5
    shell:
        "bwa-mem2 mem -M -t {threads} -R {params.rg} {input.ref} {input.r} > {output}"
            # -M = mark short split hits as secondary (for Picard compatibility)
            # -t = number of threads
            # -R = Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM.
            #      The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’.read group
            #      header line suck as '@RG\tID:foo\tSM:bar'

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 4 - Altera formato de .sam para .bam
rule view:
        input: os.path.join(config['align_path'], "{ids}.sam")
        output: temp(os.path.join(config['align_path'], "{ids}.bam"))
        wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
        threads: 5
        shell: "samtools view -b -@{threads} {input} > {output}"
        # Views and converts SAM/BAM/CRAM files.
            # -b (--bam) = Output in the BAM format.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 5 - Organiza o .bam por nome dos reads
rule sort_by_name:
        input: os.path.join(config['align_path'], "{ids}.bam")
        output: temp(os.path.join(config['align_path'], "{ids}.nameSrt.bam"))
        wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
        threads: 5
        shell: "samtools sort -n -@{threads} {input} > {output}"
        # Sort alignments by leftmost coordinates, or by read name when -n is used.
            # -n = Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
            # -@ = Set number of sorting and compression threads. By default, operation is single-threaded.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 6 - Fill in mate coordinates, ISIZE and mate related flags from a name-sorted or name-collated alignment
rule fixmate:
    input: os.path.join(config['align_path'], "{ids}.nameSrt.bam")
    output: temp(os.path.join(config['align_path'], "{ids}.fm.bam"))
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
    threads: 5
    shell:
        "samtools fixmate -m -c -@{threads} -O bam,level=1 {input} {output}"
        # Fill in mate coordinates, ISIZE and mate related flags from a name-sorted or name-collated alignment.
            # -O bam = output format
            # level=1 ?
            # -m = add mate score tags. These are used by markdup to select the best reads to keep.
            # -c = add template cigar ct tag.

#-------------------------------------------------------------------------------------------------------------------------------------------------
def join_input(wildcards):
    inputs = grouped_fqs[wildcards.ids]
    return [config['align_path']+i+'.fm.bam' for i in inputs]

# 7 - Junta os .fm.bam de um mesmo grupo em apenas um arquivo joined.bam. 
#     Exemplo: ['18_190808_L001.fm.bam', '18_190925_L001.fm.bam', '18_190930_L001.fm.bam'] serao agrupados em um arquivo 18_L001_joined.bam
rule join_bams:
    input: join_input
    output: temp(os.path.join(config['align_path'], "{ids}_joined.bam"))
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + list(grouped_fqs.keys()))
    shell:
        "samtools merge {output} {input}"
        # Merge sorted alignment files, producing a single sorted output file that contains all the input records and maintains the existig sort order.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 8 - Ordena os alignments por coordenada
rule sortbam:
    input: os.path.join(config['align_path'], "{ids}_joined.bam")
    output: temp(os.path.join(config['align_path'], "{ids}.srt.bam"))
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + list(grouped_fqs.keys()))
    threads: 5
    shell:
        "samtools sort -l 1 -@{threads} -o {output} {input}"
        # Sort alignments by leftmost coordinates, or by read name when -n is used.
            # -l = Set the desired compression level for the final output file, ranging from 0 (uncompressed) or 1
            #      (fastest but minimal compression) to 9 (best compression but slowest to write), similarly to
            #      gzip(1)'s compression level setting. If -l is not used, the default compression level will apply.
            # -@ = Set number of sorting and compression threads. By default, operation is single-threaded.
            # -o = Write the final sorted output to FILE, rather than to standard output.
            # -T = Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an existing directory, to
            #      PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique to this invocation of the sort command.
            #      By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or
            #      if output is to standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 9 - Mark duplicate alignments from a coordinate sorted file that has been run through samtools fixmate with the
#     -m option. This program relies on the MC and ms tags that fixmate provides.
rule markdup:
    input: os.path.join(config['align_path'], "{ids}.srt.bam")
    output: os.path.join(config['align_path'], "{ids}.md.srt.bam")
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + list(grouped_fqs.keys()))
    threads: 5
    shell:
        "samtools markdup -@{threads} {input} {output}"

            # -@ = Number of input/output compression threads to use in addition to main thread [0]. 

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 10 - Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access.
rule index_markdup:
    input: os.path.join(config['align_path'], "{ids}.md.srt.bam")
    output: os.path.join(config['align_path'], "{ids}.md.srt.bam.bai")
    wildcard_constraints: ids = "|".join(list(grouped_fqs.keys()))
    threads: 4
    shell:
        "samtools index -@{threads} {input}"
