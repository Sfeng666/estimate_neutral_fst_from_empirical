import optparse

def calc_minor_af_from_ct(count_table_pop1, count_table_pop2, out_afs_pop1, out_afs_pop2):
    # # for test purpose
    # count_table_pop1 = "../data/ms_simulation_chr2R_5_pop_ZI.ct"
    # count_table_pop2 = "../data/ms_simulation_chr2R_5_pop_FR.ct"
    # out_fst = "../data/fst_chr2R.txt"

    # Read the genotype count table and calculate SNP Fst for each site
    with open(count_table_pop1, 'r') as f_ct_pop1, open(count_table_pop2, 'r') as f_ct_pop2, open(out_afs_pop1, 'w') as f_out_afs_pop1, open(out_afs_pop2, 'w') as f_out_afs_pop2:
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
                # confirm the minor allele as the allele with the second highest total allele frequency across populations
                combined_af = {idx: sum([p1_afs[idx], p2_afs[idx]]) for idx in range(len(line_ct_pop1))}
                minor_af_idx = sorted(combined_af, key=combined_af.get, reverse = True)[1]
                minor_af_pop1 = min(p1_afs[minor_af_idx], 1 - p1_afs[minor_af_idx])  # since the minor allele across populations may not be the minor allele of each population, we take the smaller one of top two allele frequencies as MAF per population
                minor_af_pop2 = min(p2_afs[minor_af_idx], 1 - p2_afs[minor_af_idx])

                f_out_afs_pop1.write(f"{minor_af_pop1}\n")
                f_out_afs_pop2.write(f"{minor_af_pop2}\n")  

def main():
    usage = "usage: %prog [options] args"
    description = '''Calculate SNP Fst (Reynolds' modification to Wright's original formulation) between two populations'''
    version = '%prog 12.27.2024'
    parser = optparse.OptionParser(usage=usage,version=version, description = description)
    parser.add_option('--in_ct_pop1', dest='count_table_pop1', help='Input count table file for population 1 (.txt)', metavar = "PATH")
    parser.add_option('--in_ct_pop2', dest='count_table_pop2', help='Input count table file for population 2 (.txt)', metavar = "PATH")
    parser.add_option('--out_afs_pop1', dest='out_afs_pop1', help='Output file of per-site minor allele frequency in population 1, where each line match with the count table', metavar = "PATH")
    parser.add_option('--out_afs_pop2', dest='out_afs_pop2', help='Output file of per-site minor allele frequency in population 2, where each line match with the count table', metavar = "PATH")

    (options, args) = parser.parse_args()

    if not all([options.count_table_pop1, options.count_table_pop2, options.out_afs_pop1, options.out_afs_pop2]):
        parser.error('Missing required arguments')

    calc_minor_af_from_ct(options.count_table_pop1, options.count_table_pop2, options.out_afs_pop1, options.out_afs_pop2)

if __name__ == '__main__':
    main()