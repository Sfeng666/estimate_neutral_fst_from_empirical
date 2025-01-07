import os
import random

# Define input and output directories
pops = ["ZI", "FR"]
sample_size = 18
input_dir = "../data/input"
output_dir = "../data/output"
bed_rmrecomb_syn = os.path.join(output_dir, "all_sites_chr_dm6_rmlowrecomb_syn.bed")
in_ct_table = os.path.join(input_dir, "AllCounts_{0}_NoInv_{1}.txt")  # input count table (per chromosome) file name template
out_ct_table = os.path.join(output_dir, "AllCounts_{0}_NoInv_{1}_rmlowrecomb_syn.txt")    # output count table (per chromosome) file name template
combined_output = os.path.join(output_dir, "AllCounts_{}_NoInv_combined_rmlowrecomb_syn.txt")   # combined output count table (across chromosomes) file name template

# Create output directory if it doesn't existf
os.makedirs(output_dir, exist_ok=True)

# Function to extract 4-fold synonymous sites for a given chromosome
def extract_sites(chr):
    in_ct_table_pop1 = in_ct_table.format(pops[0], chr.replace('chr', 'Chr'))   # set file path to the original count table for pop 1
    in_ct_table_pop2 = in_ct_table.format(pops[1], chr.replace('chr', 'Chr'))   # set file path to the original count table for pop 2
    out_ct_table_pop1 = out_ct_table.format(pops[0], chr)   # set file path to the output count table for pop 1
    out_ct_table_pop2 = out_ct_table.format(pops[1], chr)   # set file path to the output count table for pop 2

    # extract rows of count table whose genomic positions match those of 4-fold synonymous sites with low recombination rates
    i = 0
    with open(in_ct_table_pop1, 'r') as f_in_ct_table_pop1, open(in_ct_table_pop2, 'r') as f_in_ct_table_pop2, open(out_ct_table_pop1, 'w') as f_out_ct_table_pop1, open(out_ct_table_pop2, 'w') as f_out_ct_table_pop2:
        for line1, line2 in zip(f_in_ct_table_pop1, f_in_ct_table_pop2):
            if all([line1, line2]):
                # try:
                #     genomic_position = pos_ct_tab[chr][i]
                # except IndexError:
                #     print(f"Reached end of chromosome {chr}. The length of chromosome {chr} is {len(pos_ct_tab[chr])}, but the count table has {i} rows.")
                #     break
                if i in pos_rmrecomb_syn[chr]:
                    # also filter for sites that have a total allele count >= sample_size
                    allele_counts_pop1 = list(int(x) for x in line1.strip().split())
                    allele_counts_pop2 = list(int(x) for x in line2.strip().split())
                    if sum(allele_counts_pop1) >= sample_size and sum(allele_counts_pop2) >= sample_size:   # only keep sites with total allele count >= sample_size in both populations for Fst calculation
                        # if allele count exceeds sample size, downsample the count for each allele proportionally without replacement, so that the total count = sample_size
                        # for pop1
                        if sum(allele_counts_pop1) == sample_size:
                            f_out_ct_table_pop1.write(line1)
                        elif sum(allele_counts_pop1) > sample_size:
                            alleles = []
                            for count, allele in zip(allele_counts_pop1, range(len(allele_counts_pop1))):
                                alleles.extend([allele] * count)
                            downsampled_alleles = random.sample(alleles, k=sample_size)
                            downsampled_counts = [downsampled_alleles.count(allele) for allele in range(len(allele_counts_pop1))]
                            f_out_ct_table_pop1.write('\t'.join(str(x) for x in downsampled_counts) + '\n')
                        # for pop2
                        if sum(allele_counts_pop2) == sample_size:
                            f_out_ct_table_pop2.write(line2)
                        elif sum(allele_counts_pop2) > sample_size:
                            alleles = []
                            for count, allele in zip(allele_counts_pop2, range(len(allele_counts_pop2))):
                                alleles.extend([allele] * count)
                            downsampled_alleles = random.sample(alleles, k=sample_size)
                            downsampled_counts = [downsampled_alleles.count(allele) for allele in range(len(allele_counts_pop2))]
                            f_out_ct_table_pop2.write('\t'.join(str(x) for x in downsampled_counts) + '\n')
            i += 1

    # Add the per-chromosome output to a combined output for each population
    os.system(f"cat {out_ct_table_pop1} >> {combined_output.format(pops[0])}")
    os.system(f"cat {out_ct_table_pop2} >> {combined_output.format(pops[1])}")

# Read positions of 4-fold synonymous sites with low recombination rates into a dictionary
chromosomes = set()
pos_rmrecomb_syn = {}
with open(bed_rmrecomb_syn, 'r') as f:
    for line in f:
        fields = line.strip().split()
        chr = fields[0]
        pos_ct = int(fields[3])     # the relative position of the site in the count table (equivalent to the genomic coordinate in the older assembly)
        if chr not in pos_rmrecomb_syn:
            pos_rmrecomb_syn[chr] = set()
        pos_rmrecomb_syn[chr].add(pos_ct)
        chromosomes.add(chr)

for chr in chromosomes:
    extract_sites(chr)

print(f"Extraction complete. Output files are located in {output_dir}.")