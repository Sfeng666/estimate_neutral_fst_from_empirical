import optparse

# # test purpose
# p1_afs = {0: 17/18, 1: 1/18}
# p2_afs = {0: 18/18, 1: 0/18}
# p1_afs = {0: 5/18, 1: 13/18}
# p2_afs = {0: 3/18, 1: 15/18}
# size1 = 18
# size2 = 18

## function to calculate SNP Reynolds Fst
def Fst_reynolds(p1_afs, p2_afs, size1, size2):

    # p1_afs and p2_afs are dicitonaries that include allele frequency spectrum in each population
    # size 1 and size 2 are allelic sample size of each population (the number of auto/X chromosomes, instead of the number of individuals). This is to simplify calculation, especially for the cases of X chromosome when a mixed gender of individuals are sampled.
    # simply implementing the Reynolds's formulation (Reynolds et al, 1983), while the variable names is indicating the location of terms in the formulations.

    # convert the input allelic (haploid) sample size to a individual (diploid) sample size
    size1 = size1/2
    size2 = size2/2

    # the first term that is shared by both numerator (al) and denominator (al + bl)
    first_term = sum(list((p1_afs[x] - p2_afs[x])**2 for x in p1_afs))/2

    # the expression at the numerator of the second term that is also shared by both numerator (al) and denominator (al + bl)
    second_term_num2 = size1*(1 - sum(list(x**2 for x in p1_afs.values()))) + size2*(1 - sum(list(x**2 for x in p2_afs.values())))

    # implement al
    al_frac_num = ((size1 + size2))*second_term_num2
    al_frac_denom = 4*size1*size2*(size1 + size2- 1)
    frac = al_frac_num/al_frac_denom
    al = first_term - frac
    # if al < 0:    # FST is allowed to be slightly negative
    #     al = 0
 
    # implement al + bl
    albl_frac_num = (4*size1*size2 - (size1 + size2))*second_term_num2
    albl_frac = albl_frac_num/al_frac_denom
    albl = first_term + albl_frac

    # calculate Fst as a ratio of al and al + bl
    try:
        snp_fst = al/albl
    except:
        snp_fst = 0 # for sites without variation in empirical data, the denominator is 0, so the Fst is set to 0
        print(p1_afs, p2_afs, al, albl)
    return snp_fst

    # # seperately return al and al + bl for Fst calculation (weighted average of sites within windows)
    # return al, albl

def calc_snp_fst_reynolds_from_ct(count_table_pop1, count_table_pop2, out_fst):
    # # for test purpose
    # count_table_pop1 = "../data/ms_simulation_chr2R_5_pop_ZI.ct"
    # count_table_pop2 = "../data/ms_simulation_chr2R_5_pop_FR.ct"
    # out_fst = "../data/fst_chr2R.txt"

    # Read the genotype count table and calculate SNP Fst for each site
    with open(count_table_pop1, 'r') as f_ct_pop1, open(count_table_pop2, 'r') as f_ct_pop2, open(out_fst, 'w') as f_out_fst:
        lines_ct_pop1 = f_ct_pop1.readlines()
        lines_ct_pop2 = f_ct_pop2.readlines()

        size1 = sum(list(int(ct) for ct in lines_ct_pop1[0].strip().split("\t")))
        size2 = sum(list(int(ct) for ct in lines_ct_pop2[0].strip().split("\t")))
        for line_ct_pop1, line_ct_pop2 in zip(lines_ct_pop1, lines_ct_pop2):
            line_ct_pop1 = line_ct_pop1.strip().split("\t")
            line_ct_pop2 = line_ct_pop2.strip().split("\t")
            p1_afs = {idx: float(line_ct_pop1[idx])/size1 for idx in range(len(line_ct_pop1))}
            p2_afs = {idx: float(line_ct_pop2[idx])/size2 for idx in range(len(line_ct_pop2))}
            if not all([(1 - sum(list(x**2 for x in p1_afs.values()))) == 0, (1 - sum(list(x**2 for x in p2_afs.values()))) == 0, p1_afs == p2_afs]):
            # if p1_afs != p2_afs:    # skip sites without variation
                snp_fst = Fst_reynolds(p1_afs, p2_afs, size1, size2)
                f_out_fst.write(f"{snp_fst}\n")

def main():
    usage = "usage: %prog [options] args"
    description = '''Calculate SNP Fst (Reynolds' modification to Wright's original formulation) between two populations'''
    version = '%prog 12.27.2024'
    parser = optparse.OptionParser(usage=usage,version=version, description = description)
    parser.add_option('--in_ct_pop1', dest='count_table_pop1', help='Input count table file for population 1 (.txt)', metavar = "PATH")
    parser.add_option('--in_ct_pop2', dest='count_table_pop2', help='Input count table file for population 2 (.txt)', metavar = "PATH")
    parser.add_option('--out_fst', dest='out_fst', help='Output file of per-site FST, where each line match with the count table', metavar = "PATH")

    (options, args) = parser.parse_args()

    if not all([options.count_table_pop1, options.count_table_pop2, options.out_fst]):
        parser.error('Missing required arguments')

    calc_snp_fst_reynolds_from_ct(options.count_table_pop1, options.count_table_pop2, options.out_fst)

if __name__ == '__main__':
    main()