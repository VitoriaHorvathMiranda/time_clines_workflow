# CORRIGIR config['raw_fqs_path'] e align_path COM O CAMINHO CERTO PARA AS AMOSTRAS FRICANAS E EUROPEIAS 
#-----------------------------------------------------------------------------------------------------------------
rule download_files_anc: #download ncbi sequences based on SRA access number (see config file)
    input: 
    output: temp(os.path.join(config['anc_folder'], "raw/{SRAs}.sra"))
    #params: out_path = os.path.join(config['out_dir_test'], "{SRAs}")
    wildcard_constraints: SRAs="|".join(list(config['haploid_SRAs_dict'].values()))
    group: "get_anc"
    log: "logs/prefetch_{SRAs}.log"
    shell: "prefetch -p -o {output} {wildcards.SRAs}"
#-----------------------------------------------------------------------------------------------------------------
rule transform_files_anc: #download ncbi sequences based on SRA access number (see config file)
    input: os.path.join(config['anc_folder'], "raw/{SRAs}.sra")
    output: r1=temp(os.path.join(config['anc_folder'], "raw/{SRAs}_1.fastq.gz")), r2=temp(os.path.join(config['anc_folder'], "raw/{SRAs}_2.fastq.gz"))
    params: out_path = os.path.join(config['anc_folder'], "raw/")
    wildcard_constraints: SRAs="|".join(list(config['haploid_SRAs_dict'].values())), r=["1","2"]
    group: "get_anc"
    log: "logs/fasterq_dump_{SRAs}.log"
    #threads: 10
    shell: "fastq-dump --skip-technical --gzip --split-files {input} --outdir {params.out_path}"
#-----------------------------------------------------------------------------------------------------------------
rule rename_files_anc:
    wildcard_constraints: sample='|'.join([re.escape(x) for x in config['haploid_SRAs_dict'].keys()]), r=["1","2"]
    input: lambda wcs: os.path.join(config['anc_folder'], "raw/%s_{r}.fastq.gz") % config['haploid_SRAs_dict'][wcs.sample]
    output:
        temp(os.path.join(config['anc_folder'], "raw/{sample}_R{r}.fastq.gz"))
    group: "get_anc"
    shell:
        "mv {input} {output}"

#-----------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 1 - Trim/Elimina reads de baixa qualidade, e junta os paired-end read (R1.fastq.gz e R2.fastq.gz) em um arquivo soh merged.fastq.gz
rule fastp_anc:
    input: r1=os.path.join(config['anc_folder'], "raw/{ids}_R1.fastq.gz"), r2=os.path.join(config['anc_folder'], "raw/{ids}_R2.fastq.gz")
    output: merged_fq=temp(os.path.join(config['anc_folder'], "raw/{ids}_merged.fastq.gz")), html=os.path.join(config['anc_folder'], "raw/{ids}_merged.fastq.gz.html")
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 4
    group: "get_anc"
    shell:
        "fastp --thread {threads} -i {input.r1} -I {input.r2} -m --include_unmerged --merged_out {output.merged_fq} -h {output.html} -R \"{wildcards.ids}\""

#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
rule ungz_anc:
    input: r=os.path.join(config['anc_folder'], "raw/{ids}_merged.fastq.gz")
    output: temp(os.path.join(config['anc_folder'], "raw/{ids}_merged.fastq"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 4
    group: "get_anc"
    shell:
        "pigz -f -p {threads} -d {input} > {output}"
        # Compress or expand files
            # -f (--force) = Force overwrite, compress .gz, links, and to terminal.
            # -p (--processes) = Allow up to n processes (default is the number of online processors).
            # -d (--decompress or --uncompress) = Decompress the compressed input.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 3.1 - Mapeia os reads contra a referencia - cria os arquivos .sam
rule bwa_anc:
    input: r=os.path.join(config['anc_folder'], "raw/{ids}_merged.fastq"), refi=config['ref_index_path'], ref=config['ref_path']
    #params: rg = lambda wildcards: read_groups[wildcards.ids]
    output: temp(os.path.join(config['anc_folder'], "align/{ids}.sam"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell:
        "bwa-mem2 mem -M -t {threads} -R \"@RG\\tID:{wildcards.ids}\\tSM:{wildcards.ids}\\tPL:Illumina\"  {input.ref} {input.r} > {output}"
            
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 4 - Altera formato de .sam para .bam
rule view_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.sam")
    output: temp(os.path.join(config['anc_folder'], "align/{ids}.bam"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell: "samtools view -b -@{threads} {input} > {output}"
        
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 5 - Organiza o .bam por nome dos reads
rule sort_by_name_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.bam")
    output: temp(os.path.join(config['anc_folder'], "align/{ids}.nameSrt.bam"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell: "samtools sort -n -@{threads} {input} > {output}"
    
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 6 - Fill in mate coordinates, ISIZE and mate related flags from a name-sorted or name-collated alignment
rule fixmate_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.nameSrt.bam")
    output: temp(os.path.join(config['anc_folder'], "align/{ids}.fm.bam"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell:
        "samtools fixmate -m -c -@{threads} -O bam,level=1 {input} {output}"
    
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 8 - Ordena os alignments por coordenada
rule sortbam_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.fm.bam")
    output: temp(os.path.join(config['anc_folder'], "align/{ids}.srt.bam"))
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell:
        "samtools sort -l 1 -@{threads} -o {output} {input}"
 
#-------------------------------------------------------------------------------------------------------------------------------------------------
# 9 - Mark duplicate alignments from a coordinate sorted file that has been run through samtools fixmate with the
#     -m option. This program relies on the MC and ms tags that fixmate provides.
rule markdup_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.srt.bam")
    output: os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam")
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 5
    group: "get_anc"
    shell:
        "samtools markdup -@{threads} {input} {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 10 - Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access.
rule index_markdup_anc:
    input: os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam")
    output: os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam.bai")
    wildcard_constraints:  ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    threads: 4
    group: "get_anc"
    shell:
        "samtools index -@{threads} {input}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule GVCFs:
    input: 
        bams = os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam"), 
        bais = os.path.join(config['anc_folder'], "align/{ids}.md.srt.bam.bai"),
        ref = config['ref_path'],
    output:  os.path.join(config['anc_folder'], "call/{ids}_gatk.g.vcf.gz")
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    shell:
        "gatk --java-options '-Xmx4g' HaplotypeCaller -R {input.ref} -I {input.bams} --min-base-quality-score 20 --minimum-mapping-quality 20 --sample-ploidy 1 -O {output} -ERC GVCF"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule make_map_file:
    input: expand(os.path.join(config['anc_folder'], "call/{ids}_gatk.g.vcf.gz"), ids = config['haploid_SRAs_dict'].keys())
    output: os.path.join(config['anc_folder'], "call/gvcfs.sample_map")
    wildcard_constraints: ids = "|".join(list(config['haploid_SRAs_dict'].keys()))
    shell:
        "paste <(echo {input} | sed 's/ /\\n/g' | sort | cut -d '/' -f 7 | cut -d '_' -f 1) <(echo {input} | sed 's/ /\\n/g' | sort) > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule Genomic_DB: #não consegui rodar com o snakemake, o gatk reclama do caminho por algum motivo. Rodei direto no terminal
    input: 
        map_file = os.path.join(config['anc_folder'], "call/gvcfs.sample_map"),
        gvcfs = expand(os.path.join(config['anc_folder'], "call/{ids}_gatk.g.vcf.gz"), ids = config['haploid_SRAs_dict'].keys())
    output: os.path.join(config['anc_folder'], "call/databases/{chrom}/haploid_chrom_{chrom}_database/vcfheader.vcf")
    params: database_path = os.path.join(config['anc_folder'], "call/databases/{chrom}/haploid_chrom_{chrom}_database") # --genomicsdb-workspace-path must point to a non-existent or empty directory
    wildcard_constraints: 
        ids = "|".join(list(config['haploid_SRAs_dict'].keys())),
        chrom = "|".join(config['chrom'])
    shell: "gatk --java-options '-Xmx4g -Xms4g' GenomicsDBImport --genomicsdb-workspace-path {params.database_path} --sample-name-map {input.map_file} -L {wildcards.chrom}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
rule Genotype_call: #rodei direto no terminal, mas acho que arrumei o problema por aqui
    input: 
        database = os.path.join(config['anc_folder'], "call/databases/{chrom}/haploid_chrom_{chrom}_database/vcfheader.vcf"),
        ref = config['ref_path']
    output: os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.vcf.gz")
    params: database_path = "gendb://" + os.path.join(config['anc_folder'], "call/databases/{chrom}/haploid_chrom_{chrom}_database")
    wildcard_constraints: chrom = "|".join(config['chrom'])
    shell: "gatk --java-options '-Xmx4g' GenotypeGVCFs -R {input.ref} -V {params.database_path} -O {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Remove os indels e os SNPs multialélicos
rule remove_indels:
    input: 
        ref = config['ref_path'],
        vcf = os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.vcf.gz")
    output: temp(os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.biallelic.vcf.gz"))
    shell: "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------
# filter repeat regions from vcf
rule filter_repeat_anc:
    input: 
        vcf = os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.biallelic.vcf.gz"),
        bed_file = os.path.join(config['ref_folder'],"repeat_6.bed")
    output: os.path.join(config['anc_folder'], "call/AFR_EU_chrom_{chrom}.biallelic.clean.vcf")
    shell: "bedtools subtract -header -a {input.vcf} -b {input.bed_file} > {output}"

#-------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------------
#rule make_sync_anc: #FALTA ARRUMAR
#    input:
#        ref=config['ref_path'],
#        merged_sync = os.path.join(config['analysis_path'], "ancestry/all_samples.sync"),
#        bam_path = config['']
#    output:
#    params: bam_path = config[''], outdir = config['']
#    shell: "grenedalf sync \
#    --sam-path {params.bam_path} \
#    --sam-min-map-qual 20 \
#    --sam-min-base-qual 20 \
#    --reference-genome-fasta {input.ref} \
#    --filter-region-vcf {input.merged_sync} \
#    --out-dir {params.outdir}"