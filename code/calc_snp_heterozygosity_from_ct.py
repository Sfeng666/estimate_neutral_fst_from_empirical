import optparse

def calc_total_heterozygosity_from_ct(count_table_pop1, count_table_pop2, out_heterozygosity):
    # # for test purpose
    # count_table_pop1 = "../data/ms_simulation_chr2R_5_pop_ZI.ct"
    # count_table_pop2 = "../data/ms_simulation_chr2R_5_pop_FR.ct"
    # out_fst = "../data/fst_chr2R.txt"

    # Read the genotype count table and calculate SNP Fst for each site
    with open(count_table_pop1, 'r') as f_ct_pop1, open(count_table_pop2, 'r') as f_ct_pop2, open(out_heterozygosity, 'w') as f_out_heterozygosity:
        lines_ct_pop1 = f_ct_pop1.readlines()
        lines_ct_pop2 = f_ct_pop2.readlines()

        size1 = sum(list(int(ct) for ct in lines_ct_pop1[0].strip().split("\t")))
        size2 = sum(list(int(ct) for ct in lines_ct_pop2[0].strip().split("\t")))
        for line_ct_pop1, line_ct_pop2 in zip(lines_ct_pop1, lines_ct_pop2):
            line_ct_pop1 = line_ct_pop1.strip().split("\t")
            line_ct_pop2 = line_ct_pop2.strip().split("\t")
            p1_afs = {idx: float(line_ct_pop1[idx])/size1 for idx in range(len(line_ct_pop1))}
            p2_afs = {idx: float(line_ct_pop2[idx])/size2 for idx in range(len(line_ct_pop2))}
            if not all([(1 - sum(list(x**2 for x in p1_afs.values()))) == 0, (1 - sum(list(x**2 for x in p2_afs.values()))) == 0, p1_afs == p2_afs]):   # skip sites that are fixed at the same allele in both populations
                total_heterozygosity = (1 - sum(list(x**2 for x in p1_afs.values()))) + (1 - sum(list(x**2 for x in p2_afs.values())))

                f_out_heterozygosity.write(f"{total_heterozygosity}\n")

def main():
    usage = "usage: %prog [options] args"
    description = '''Calculate SNP Fst (Reynolds' modification to Wright's original formulation) between two populations'''
    version = '%prog 12.27.2024'
    parser = optparse.OptionParser(usage=usage,version=version, description = description)
    parser.add_option('--in_ct_pop1', dest='count_table_pop1', help='Input count table file for population 1 (.txt)', metavar = "PATH")
    parser.add_option('--in_ct_pop2', dest='count_table_pop2', help='Input count table file for population 2 (.txt)', metavar = "PATH")
    parser.add_option('--out_heterozygosity', dest='out_heterozygosity', help='Output file of per-site total heterozygosity across populations, where each line match with the count table', metavar = "PATH")

    (options, args) = parser.parse_args()

    if not all([options.count_table_pop1, options.count_table_pop2, options.out_heterozygosity]):
        parser.error('Missing required arguments')

    calc_total_heterozygosity_from_ct(options.count_table_pop1, options.count_table_pop2, options.out_heterozygosity)

if __name__ == '__main__':
    main()