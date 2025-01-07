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

    # # calculate Fst as a ratio of al and al + bl
    # try:
    #     snp_fst = al/albl
    # except:
    #     print(p1_afs, p2_afs)
    # return snp_fst

    # seperately return al and al + bl for Fst calculation (weighted average of sites within windows)
    return al, albl

## function to calculate per-site heterozygosity
def calc_heterozygosity(afs):
    heterozygosity = 1 - sum(list(afs[x]**2 for x in afs))
    return heterozygosity

def calc_win_fst_diversity_from_ct(win_size, win_info, count_table_pop1, count_table_pop2, out_fst, out_diversity_pop1, out_diversity_pop2):
    # # for test purpose
    # count_table_pop1 = "../data/ms_simulation_chr2R_5_pop_ZI.ct"
    # count_table_pop2 = "../data/ms_simulation_chr2R_5_pop_FR.ct"
    # win_info = "../data/ms_simulation_chr2R_5_pop.win"
    # out_fst = "../data/fst_chr2R_5_pop.txt"

    # Read the genotype count table and calculate SNP Fst for each site
    with open(count_table_pop1, 'r') as f_ct_pop1, open(count_table_pop2, 'r') as f_ct_pop2, open(win_info, 'r') as f_win, open(out_fst, 'w') as f_out_fst, open(out_diversity_pop1, 'w') as f_out_diversity_pop1, open(out_diversity_pop2, 'w') as f_out_diversity_pop2:
        lines_ct_pop1 = f_ct_pop1.readlines()
        lines_ct_pop2 = f_ct_pop2.readlines()
        num_segsites_allwins = list(int(line.strip()) for line in f_win.readlines())    # list of the number of segregating sites within each window

        size1 = sum(list(int(ct) for ct in lines_ct_pop1[0].strip().split("\t")))
        size2 = sum(list(int(ct) for ct in lines_ct_pop2[0].strip().split("\t")))

        current_win = 0  # index of the current window
        num_segsites_currentwin = num_segsites_allwins[current_win]   # number of segregating sites within the first window
        num_segsites_calculated = 0
        al_currentwin = []  # list of al of all sites within the current window
        albl_currentwin = []    # list of al + bl of all sites within the current window
        p1_heterozygosity_currentwin = []  # list of heterozygosity of all sites within the current window for population 1
        p2_heterozygosity_currentwin = []  # list of heterozygosity of all sites within the current window for population 2

        for line_ct_pop1, line_ct_pop2 in zip(lines_ct_pop1, lines_ct_pop2):
            line_ct_pop1 = line_ct_pop1.strip().split("\t")
            line_ct_pop2 = line_ct_pop2.strip().split("\t")
            p1_afs = {idx: float(line_ct_pop1[idx])/size1 for idx in range(len(line_ct_pop1))}
            p2_afs = {idx: float(line_ct_pop2[idx])/size2 for idx in range(len(line_ct_pop2))}

            # for FST
            al, albl = Fst_reynolds(p1_afs, p2_afs, size1, size2)   # calculate the numerator and denominator for per-site FST
            al_currentwin.append(al)
            albl_currentwin.append(albl)

            # for heterozygosity
            p1_heterozygosity = calc_heterozygosity(p1_afs)
            p2_heterozygosity = calc_heterozygosity(p2_afs)
            p1_heterozygosity_currentwin.append(p1_heterozygosity)
            p2_heterozygosity_currentwin.append(p2_heterozygosity)

            num_segsites_calculated += 1
            if num_segsites_calculated == num_segsites_currentwin:
                # calculate window FST as the weighted average of al and al + bl within the current window
                fst_currentwin = sum(al_currentwin)/sum(albl_currentwin)
                f_out_fst.write(f"{fst_currentwin}\n")
                al_currentwin = []  # reset the list of al for the next window
                albl_currentwin = []    # reset the list of al + bl for the next window

                # calculate nucleotide diversity of the current window for each population
                diversity_currentwin_pop1 = sum(p1_heterozygosity_currentwin)/int(win_size)
                diversity_currentwin_pop2 = sum(p2_heterozygosity_currentwin)/int(win_size)
                f_out_diversity_pop1.write(f"{diversity_currentwin_pop1}\n")
                f_out_diversity_pop2.write(f"{diversity_currentwin_pop2}\n")
                p1_heterozygosity_currentwin = []  # reset the list of heterozygosity for the next window
                p2_heterozygosity_currentwin = []  # reset the list of heterozygosity for the next window

                # move to the next window, if there is any
                if current_win < len(num_segsites_allwins) - 1:
                    current_win += 1
                    num_segsites_currentwin = num_segsites_allwins[current_win]
                    num_segsites_calculated = 0

def main():
    usage = "usage: %prog [options] args"
    description = '''Calculate window Fst (Reynolds' estimator) and nucleotide diversity between two populations'''
    version = '%prog 12.20.2024'
    parser = optparse.OptionParser(usage=usage,version=version, description = description)
    parser.add_option('--win_size', dest='win_size', help='Input parameter for window size (bp)', metavar = "INT")
    parser.add_option('--in_win', dest='win_info', help='Input table of the number of segregating sites across populations within each window (.win)', metavar = "PATH")
    parser.add_option('--in_ct_pop1', dest='count_table_pop1', help='Input count table file for population 1 (.ct)', metavar = "PATH")
    parser.add_option('--in_ct_pop2', dest='count_table_pop2', help='Input count table file for population 2 (.ct)', metavar = "PATH")
    parser.add_option('--out_fst', dest='out_fst', help='Output file of per-site FST, where each line match with the count table', metavar = "PATH")
    parser.add_option('--out_diversity_pop1', dest='out_diversity_pop1', help='Output file of window nucleotide diversity for population 1', metavar = "PATH")
    parser.add_option('--out_diversity_pop2', dest='out_diversity_pop2', help='Output file of window nucleotide diversity for population 2', metavar = "PATH")

    (options, args) = parser.parse_args()

    if not all([options.win_size, options.win_info, options.count_table_pop1, options.count_table_pop2, options.out_fst, options.out_diversity_pop1, options.out_diversity_pop2]):
        parser.error('Missing required arguments')

    calc_win_fst_diversity_from_ct(options.win_size, options.win_info, options.count_table_pop1, options.count_table_pop2, options.out_fst, options.out_diversity_pop1, options.out_diversity_pop2)

if __name__ == '__main__':
    main()


# ## function to calculate SNP Reynolds Fst
# def Fst_reynolds(p1_afs, p2_afs, size1, size2):
#     first_term = sum(list((p1_afs[x] - p2_afs[x])**2 for x in p1_afs))/2
#     second_term_num2 = size1*(1 - sum(list(x**2 for x in p1_afs.values()))) + size2*(1 - sum(list(x**2 for x in p2_afs.values())))
#     al_frac_num = ((size1 + size2))*second_term_num2
#     al_frac_denom = 4*size1*size2*(size1 + size2- 1)
#     frac = al_frac_num/al_frac_denom
#     al = first_term - frac
#     albl_frac_num = (4*size1*size2 - (size1 + size2))*second_term_num2
#     albl_frac = albl_frac_num/al_frac_denom
#     albl = first_term + albl_frac
#     return al, albl

# # for test purpose
# count_table_pop1 = "../data/ms_simulation_chr2R_5_pop_ZI.ct"
# count_table_pop2 = "../data/ms_simulation_chr2R_5_pop_FR.ct"
# win_info = "../data/ms_simulation_chr2R_5_pop.win"
# out_fst = "../data/fst_chr2R.txt"

# # Read the genotype count table and calculate SNP Fst for each site
# with open(count_table_pop1, 'r') as f_ct_pop1, open(count_table_pop2, 'r') as f_ct_pop2, open(win_info, 'r') as f_win, open(out_fst, 'w') as f_out_fst:
#     lines_ct_pop1 = f_ct_pop1.readlines()
#     lines_ct_pop2 = f_ct_pop2.readlines()
#     num_segsites_allwins = list(int(line.strip()) for line in f_win.readlines())    # list of the number of segregating sites within each window
#     size1 = sum(list(int(ct) for ct in lines_ct_pop1[0].strip().split("\t")))
#     size2 = sum(list(int(ct) for ct in lines_ct_pop2[0].strip().split("\t")))
#     current_win = 0  # index of the current window
#     num_segsites_currentwin = num_segsites_allwins[current_win]   # number of segregating sites within the first window
#     num_segsites_calculated = 0
#     al_currentwin = []  # list of al values within the current window
#     albl_currentwin = []    # list of al + bl values within the current window
#     for line_ct_pop1, line_ct_pop2 in zip(lines_ct_pop1, lines_ct_pop2):
#         line_ct_pop1 = line_ct_pop1.strip().split("\t")
#         line_ct_pop2 = line_ct_pop2.strip().split("\t")
#         p1_afs = {idx: float(line_ct_pop1[idx])/size1 for idx in range(len(line_ct_pop1))}
#         p2_afs = {idx: float(line_ct_pop2[idx])/size2 for idx in range(len(line_ct_pop2))}
#         al, albl = Fst_reynolds(p1_afs, p2_afs, size1, size2)
#         al_currentwin.append(al)
#         albl_currentwin.append(albl)
#         num_segsites_calculated += 1
#         if num_segsites_calculated == num_segsites_currentwin:
#             # calculate the weighted average of al and al + bl within the current window
#             fst_currentwin = sum(al_currentwin)/sum(albl_currentwin)
#             f_out_fst.write(f"{fst_currentwin}\n")
#             # move to the next window
#             current_win += 1
#             num_segsites_currentwin = num_segsites_allwins[current_win]
#             num_segsites_calculated = 0
