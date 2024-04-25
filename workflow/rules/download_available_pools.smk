#########################################

#-----------------------------------------------------------------------------------------------------------------
rule download_files: #download ncbi sequences based on SRA access number (see config file)
    input: 
    output: temp(os.path.join(config['out_dir_test'], "{SRAs}.sra"))
    #params: out_path = os.path.join(config['out_dir_test'], "{SRAs}")
    wildcard_constraints: SRAs="|".join(list(config['SRAs_dict'].values()))
    #threads: 10
    shell: "prefetch -p -o {output} {wildcards.SRAs}"
#-----------------------------------------------------------------------------------------------------------------
rule transform_files: #download ncbi sequences based on SRA access number (see config file)
    input: os.path.join(config['out_dir_test'], "{SRAs}.sra")
    output: r1=os.path.join(config['out_dir_test'], "{SRAs}_1.fastq.gz"), r2=os.path.join(config['out_dir_test'], "{SRAs}_2.fastq.gz")
    params: out_path = config['out_dir_test']
    wildcard_constraints: SRAs="|".join(list(config['SRAs_dict'].values())), r=["1","2"]
    threads: 10
    shell: "fasterq-dump --skip-technical --split-files --threads {threads} --gzip {input} --outdir {params.out_path}"
#-----------------------------------------------------------------------------------------------------------------
rule rename_files:
    wildcard_constraints: sample='|'.join([re.escape(x) for x in config['SRAs_dict'].keys()]), r=["1","2"]
    input: lambda wcs: os.path.join(config['out_dir_test'], "%s_{r}.fastq.gz") % config['SRAs_dict'][wcs.sample]
    output:
        os.path.join(config['raw_fqs_path'], "{sample}_L001_R{r}.fastq.gz")
    shell:
        "mv {input} {output}"



#rule rename_downloads: #deve ter um jeito melhor de fazer isso 
#    input:
#        SRAs = [os.path.join(config['out_dir_test'], SRAs + ".fastq.gz") for SRAs in config['SRAs'].keys()],
#        SRAs_info = config['SRA_info']
#    output:
#        [os.path.join(config['out_dir_test'], samples + ".fastq.gz") for samples in config['SRAs'].values()]
#    params: out_path = config['out_dir_test'] 
#    shell:
#        "bash scripts/bash/fix_sra_names.sh {input.SRAs_info} {params.out_path}" #quando acrescentar amostras vai mudar a data dos outros arquivos também?

#-----------------------------------------------------------------------------------------------------------------
#rule fastp_dl:
#    input: os.path.join(config['out_dir_test'], "{samples}.fastq.gz") 
#    output: 
#        filtered = os.path.join(config['processed_path'], "{samples}_merged.fastq.gz"), 
#        html_report = os.path.join(config['processed_path'], "{samples}_merged.fastq.gz.html")
#    wildcard_constraints: samples = config['SRAs'].values()
#    threads: 5
#    shell: "fastp --thread {threads} -i {input} -o {output.filtered} -h {output.html_report}"