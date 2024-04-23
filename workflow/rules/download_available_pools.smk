#########################################

#-----------------------------------------------------------------------------------------------------------------
rule download_files: #download ncbi sequences based on SRA access number (see config file)
    input: 
    output: os.path.join(config['out_dir_test'], "{SRAs}.fastq.gz", SRAs=config['SRAs'])
    params: out_path = config['out_dir_test']
    shell: "fastq-dump --skip-technical --gzip {wildcards.SRAs} --outdir {params.out_path}"

#-----------------------------------------------------------------------------------------------------------------
rule rename_downloads:
    input:
        os.path.join(config['out_dir_test'], "{SRAs}.fastq.gz", SRAs=config['SRAs'])
    output:
        #lambda wildcards: 
        os.path.join(config['out_dir_test'], config['SRA_names'][wildcards.SRAs] + ".fastq.gz")
    shell:
        "mv {input} {output}"     