import snakemake
import glob
import re
import gzip
import pandas as pd

ref_index_path = os.path.join(config['ref_path'], ".pac")
r1_fqs = glob.glob(os.path.join(config['raw_fqs_path'], f"*_R1.fastq.gz")) # procura todos arquivos com final _R1.fastq.gz
fqs_pref = [str.replace(str.replace(fqname, '_R1.fastq.gz', ''), config['raw_fqs_path'], '') for fqname in r1_fqs] # guarda o prefixo dos arquivos .fastq.gz encontrados anteriormente. Exemplo: “data3/murillo/time_clines/data/seqs/raw/13_190808_L001_R1.fastq.gz" é armazenado como "13_190808_L001".
grouped_fqs = {}
for fq in fqs_pref:
    finds = re.findall(r'(^.{2,5})_(19\d+|22.+)_.*L00', fq) # # guarda partes entre () do prefixo se for do padrao "**_19****_L001". Exemplo: "13_190808_L001" -> finds =  [('13', '190808')]
    assert len(finds) <= 1 # interrompe execucao se expressao for FALSE
    if len(finds) == 0: # caso nao seja do padrao "**_19****_L001". Outros padroes sao "*_L1A1_*_L001" e "*_L002".
        assert fq not in grouped_fqs # testa se prefixo fq ja nao esta em grouped_fqs
        new_fqname = fq # se passou do assert acima, eh um novo grupo
        if 'L1A1' in fq: # se prefixo possuir "L1A1"
            finds = re.findall(r'(^.{2,3})_L1A1.+(_L00\d)', fq) # guarda partes entre () do prefixo se for do padrao "*_L1A1_*_L001". Exemplo: "09_L1A1_EKD...SXX_L001" -> finds = [('09', 'L001')]
            assert len(finds) == 1 
            new_fqname = "".join(finds[0]) # junta as duas partes de finds. Exemplo: finds = [('09', 'L001')] -> new_fqname = 09_L001
        grouped_fqs[new_fqname] = [fq] # cria novo grupo new_fqname e adiciona o prefixo fq ao grupo
    else:
        joined_name = str.replace(fq, f'_{finds[0][1]}', '') # remove parte do prefixo dos que possuem padrao "**_19****_L001". Exemplo: "13_190808_L001" -> "13_L001"
        if joined_name not in grouped_fqs: # se grupo nao existe
            grouped_fqs[joined_name] = [fq] # cria novo grupo e adiciona prefixo ao grupo
        else:
            grouped_fqs[joined_name].append(fq) # se grupo existe, adiciona prefixo ao grupo

# Ao final desse bloco, teremos os seguintes grupos + prefixos:
# {'09_L001':    ['09_L1A1_EKDL190133983-2a-5-AK3819_HWKNMDSXX_L001'], 
#  '10_L001':    ['10_L1A1_EKDL190133983-2a-AK2939-AK10154_HWKNMDSXX_L001'], 
#  '11_L001':    ['11_L1A1_EKDL190133983-2a-11-AK7152_HWKNMDSXX_L001'], 
#  '12_L001':    ['12_L1A1_EKDL190133983-2a-AK2957-AK998_HWKNMDSXX_L001'], 
#  '13_L001':    ['13_190808_L001', '13_190925_L001', '13_190930_L001', '13_191003_L001', '13_191205_L001', '13_191209_L001'], 
#  '15_L001':    ['15_190808_L001'], 
#  '16_L001':    ['16_190808_L001', '16_190925_L001', '16_190930_L001', '16_191003_L001', '16_191205_L001', '16_191209_L001'], 
#  '17_L001':    ['17_190808_L001', '17_190925_L001', '17_190930_L001', '17_191003_L001', '17_191205_L001', '17_191209_L001'], 
#  '18_L001':    ['18_190808_L001', '18_190925_L001', '18_190930_L001'], 
#  'HMA17_L002': ['HMA17_L002'], 
#  'MCT09_L002': ['MCT09_L002'], 
#  'MCT17_L002': ['MCT17_L002'], 
#  'MCT97_L002': ['MCT97_L002']}

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 1 - Trim/Elimina reads de baixa qualidade, e junta os paired-end read (R1.fastq.gz e R2.fastq.gz) em um arquivo soh merged.fastq.gz
rule fastp:
    input: r1=os.path.join(config['raw_fqs_path'], "{ids}_R1.fastq.gz"), r2=os.path.join(config['raw_fqs_path'], "{ids}_R2.fastq.gz")
    output: merged_fq=temp(os.path.join(config['processed_path'], "{ids}_merged.fastq.gz")), html=os.path.join(config['processed_path'], "{ids}_merged.fastq.gz.html")
    wildcard_constraints: ids = "|".join(list(config['SRAs_dict'].keys()) + fqs_pref)
    threads: 4
    shell:
        "fastp --thread {threads} -i {input.r1} -I {input.r2} -m --include_unmerged --merged_out {output.merged_fq} -h {output.html} -R \"{wildcards.ids}\""
        # fastp is a FASTQ data pre-processing tool. The algorithm has functions for quality control, trimming of adapters, filtering by quality, and read pruning. It also supports multi-threading.
            # -i = for single-end data, you only have to specify read1 input by -i or --in1, and specify read1 output by -o or --out1.
            # -I = for paired-end data, you should also specify read2 input by -I or --in2, and specify read2 output by -O or --out2.
            #     If you don't specify the output file names, no output files will be written, but the QC will still be done for
            #     both data before and after filtering. The output will be gzip-compressed if its file name ends with .gz
            # -m = for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be
            #      written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2.
            #      The merging mode is disabled by default.
            # --merged_out = in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=]).
            # -h (--html) = the html format report file name (string [=fastp.html]).
            # -R (--report_title) = should be quoted with ' or ", default is "fastp report" (string [=fastp report]).




