import snakemake
import glob
import re
import gzip
import pandas as pd


ref_path = "/dados/time_clines/data/reference/dmel-all-chromosome-r6.25.fasta" # gene de referencia
ref_index_path = ref_path + ".pac" # index de referencia
raw_fqs_path = "/dados/time_clines/data/seqs/raw/" # local dos arquivos fastq.gz
meta_path = "/dados/time_clines/data/meta/seq_metadata.tsv" # metadados
#meta = pd.read_csv(meta_path, sep="\t") # meta csv

r1_fqs = glob.glob(f"/dados/time_clines/data/seqs/raw/*_R1.fastq.gz") # procura todos arquivos com final _R1.fastq.gz
fqs_pref = [str.replace(str.replace(fqname, '_R1.fastq.gz', ''), '/dados/time_clines/data/seqs/raw/', '') for fqname in r1_fqs] # guarda o prefixo dos arquivos .fastq.gz encontrados anteriormente. Exemplo: “data3/murillo/time_clines/data/seqs/raw/13_190808_L001_R1.fastq.gz" é armazenado como "13_190808_L001".
grouped_fqs = {}
for fq in fqs_pref:
    finds = re.findall(r'(^.{2,5})_(19\d+)_.*L00', fq) # # guarda partes entre () do prefixo se for do padrao "**_19****_L001". Exemplo: "13_190808_L001" -> finds =  [('13', '190808')]
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

test=",".join(list(grouped_fqs.keys()))
#print(test)

male_ids = ["dlSC10", "dlGA10", "dlFL10", "HFL97downto60mi"]

IDs = list(grouped_fqs.keys())
IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001", "A41_L002")] #deletes samples with low quality from the following analysis
male_labels = [id + "_L001" for id in male_ids]
IDs_female_pools = [e for e in IDs_deleted_samples if e not in male_labels]


print("IDs_deleted_samples:")
print(IDs_deleted_samples)

#print("IDs_female_pools:")
#print(IDs_female_pools)
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


#IDs = list(grouped_fqs.keys())
#IDs_deleted_samples = IDs.remove("17_L001", "09_L001")
#IDs_deleted_samples = [e for e in IDs if e not in ("17_L001", "09_L001")]


#with open(meta_path, 'r') as file:
#    meta = [line.strip().split('\t') for line in file]


#list_of_ids = []
#for string in IDs_deleted_samples:
#    ids = string.split("_")[0]  # Get the part before the first "_"
#    list_of_names.append(ids)

#print(list_of_names)

#extracted_values = []
#for string in list_of_names:
#    for row in meta:
#        if string == row[0]:
#            extracted_values.append(row[1])
#            break  # Stop searching once a match is found
#    
#print(extracted_values)

