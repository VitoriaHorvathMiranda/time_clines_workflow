#########################################

#-----------------------------------------------------------------------------------------------------------------
rule download_files: #download ncbi sequences based on SRA access number (see config file)
    input: 
    output: temp(os.path.join(config['out_dir_test'], "{SRAs}.fastq.gz"))
    params: out_path = config['out_dir_test']
    wildcard_constraints: SRAs=config['SRAs'].keys()
    shell: "fastq-dump --skip-technical --gzip {wildcards.SRAs} --outdir {params.out_path}"

#-----------------------------------------------------------------------------------------------------------------
rule rename_downloads:
    input:
        SRAs = [os.path.join(config['out_dir_test'], SRAs + ".fastq.gz") for SRAs in config['SRAs'].keys()],
        SRAs_info = config['SRA_info']
    output:
        [os.path.join(config['out_dir_test'], samples + ".fastq.gz") for samples in config['SRAs'].values()]
    params: out_path = config['out_dir_test'] 
    shell:
        "bash scripts/bash/fix_sra_names.sh {input.SRAs_info} {params.out_path}"     