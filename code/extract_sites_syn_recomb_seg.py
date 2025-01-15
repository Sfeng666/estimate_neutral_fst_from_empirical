import os
import random

# Define input and output directories
pops = ["ZI", "FR"]
sample_size = 18
input_dir = "../data/input"
output_dir = "../data/output"
bed_rmrecomb_syn = os.path.join(output_dir, "all_sites_chr_dm6_rmlowrecomb_syn.bed")
in_ct_table = os.path.join(input_dir, "AllCounts_{0}_NoInv_{1}.txt")  # file name template of input count table (per chromosome) 
out_ct_table = os.path.join(output_dir, "AllCounts_{0}_NoInv_{1}_rmlowrecomb_syn_seg.txt")    # file name template of output count table of 4-fold synonymous segregating sites excluding low-recombination regions (per chromosome)
combined_output = os.path.join(output_dir, "AllCounts_{}_NoInv_combined_rmlowrecomb_syn_seg.txt")   # file name template of combined output count table (across chromosomes) 

# Create output directory if it doesn't existf
os.makedirs(output_dir, exist_ok=True)

# Overwrite 'combined_output' if it already exists
for pop in pops:
    if os.path.exists(combined_output.format(pop)):
        os.remove(combined_output.format(pop))

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
                    # also filter for sites that have a total allele count < sample_size
                    allele_counts_pop1 = list(int(x) for x in line1.strip().split())
                    allele_counts_pop2 = list(int(x) for x in line2.strip().split())
                    if sum(allele_counts_pop1) >= sample_size and sum(allele_counts_pop2) >= sample_size:   # only keep sites with total allele count >= sample_size in both populations for Fst calculation
                        # if len(set(list(idx for idx in p1_afs if p1_afs[idx] != 0)) | set(list(idx for idx in p2_afs if p2_afs[idx] != 0))) == 2:   # only keep sites with two alleles across populations
                        if not all([allele_counts_pop1.count(0) == 3, allele_counts_pop2.count(0) == 3, allele_counts_pop1 == allele_counts_pop2]):   # before downsampling, skip sites that are fixed at same allele in both populations

                            # if allele count exceeds sample size, downsample the count for each allele proportionally without replacement, so that the total count = sample_size
                            # for pop1
                            if sum(allele_counts_pop1) == sample_size:
                                out_line_pop1 = line1
                                downsampled_counts_pop1 = allele_counts_pop1
                            elif sum(allele_counts_pop1) > sample_size:
                                alleles = []
                                for count, allele in zip(allele_counts_pop1, range(len(allele_counts_pop1))):
                                    alleles.extend([allele] * count)
                                downsampled_alleles_pop1 = random.sample(alleles, k=sample_size)
                                downsampled_counts_pop1 = [downsampled_alleles_pop1.count(allele) for allele in range(len(allele_counts_pop1))]
                                out_line_pop1 = '\t'.join(str(x) for x in downsampled_counts_pop1) + '\n'
                            # for pop2
                            if sum(allele_counts_pop2) == sample_size:
                                out_line_pop2 = line2
                                downsampled_counts_pop2 = allele_counts_pop2
                            elif sum(allele_counts_pop2) > sample_size:
                                alleles = []
                                for count, allele in zip(allele_counts_pop2, range(len(allele_counts_pop2))):
                                    alleles.extend([allele] * count)
                                downsampled_alleles_pop2 = random.sample(alleles, k=sample_size)
                                downsampled_counts_pop2 = [downsampled_alleles_pop2.count(allele) for allele in range(len(allele_counts_pop2))]
                                out_line_pop2 = '\t'.join(str(x) for x in downsampled_counts_pop2) + '\n'

                            if not all([downsampled_counts_pop1.count(0) == 3, downsampled_counts_pop2.count(0) == 3, downsampled_counts_pop1 == downsampled_counts_pop2]):   # after down-sampling skip sites that are fixed at same allele in both populations
                                f_out_ct_table_pop1.write(out_line_pop1)
                                f_out_ct_table_pop2.write(out_line_pop2)
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