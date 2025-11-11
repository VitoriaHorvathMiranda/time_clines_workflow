import csv
import argparse

# Initialize the parser
parser = argparse.ArgumentParser(description="Process merged outliers and genes files and chromosomal arm")

# Add arguments
parser.add_argument("windows", type=argparse.FileType('r'), help="file with merged outlier windows")
parser.add_argument("genes", type=argparse.FileType('r'), help="file with genes")
parser.add_argument("chrom", type=str, help="chromosomal arm")
parser.add_argument("output_path", type=str, help="Path to output file")


args = parser.parse_args()

#inputchrom = "2L"
# Load genes from annot_chr##.txt
genes = []
with args.genes as gene_file:
    reader = csv.reader(gene_file, delimiter="\t")
    for row in reader:
        fbgn = row[0]
        start = int(row[4])
        end = int(row[5])
        genes.append((fbgn, start, end))

# Load regions from fst_out1percWin_regions_perChrom.tsv
regions = []
with args.windows as region_file:
    reader = csv.reader(region_file, delimiter="\t")
    #next(reader)  # Skip header
    for row in reader:
        chrom = row[0]
        if chrom == args.chrom:
            start = int(row[1])
            end = int(row[2])
            info = row[3]
            regions.append({
                'start': start,
                'end': end,
                'info': info,
                'fbgns': []
            })

# Function to check if two intervals overlap
def overlaps(a_start, a_end, b_start, b_end):
    return a_end >= b_start and b_end >= a_start

# Assign genes to regions
for fbgn, g_start, g_end in genes:
    overlapping_regions = []
    for idx, region in enumerate(regions):
        if overlaps(g_start, g_end, region['start'], region['end']):
            overlapping_regions.append(idx)

    if overlapping_regions:
        # Prefer region where gene START is within
        assigned_idx = None
        for idx in overlapping_regions:
            r = regions[idx]
            if r['start'] <= g_start <= r['end']:
                assigned_idx = idx
                break
        if assigned_idx is None:
            assigned_idx = overlapping_regions[0]  # Assign to first if no region contains the start

        regions[assigned_idx]['fbgns'].append(fbgn)

# Write result
with open(args.output_path, 'w') as out_file:
    writer = csv.writer(out_file, delimiter="\t")
    writer.writerow(["chrm","start","end","nwindows","fbgn"])
    for region in regions:
        fbgn_str = "/".join(region['fbgns'])
        writer.writerow([args.chrom, region['start'], region['end'], region['info'], fbgn_str])

