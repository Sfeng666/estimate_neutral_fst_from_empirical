#! /bin/bash
# Estimate Fst and diversity

# fixed parameters
pop1="ZI"
pop2="FR"

# demographic parameters that vary across chromosomal arms or combinations of populations
chrs=("chr2L" "chr2R" "chr3L" "chr3R" "chrX" "combined")

# summary="../data/output/estimate_heterozygosity.report"
for chr in "${chrs[@]}"; do
    in_ct_pop1="../data/output/AllCounts_${pop1}_NoInv_${chr}_rmlowrecomb_syn_seg.txt"
    in_ct_pop2="../data/output/AllCounts_${pop2}_NoInv_${chr}_rmlowrecomb_syn_seg.txt"
    out_heterozygosity="../data/output/heterozygosity_${chr}.txt"
    # estimate Fst and diversity
    python calc_snp_heterozygosity_from_ct.py --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_heterozygosity $out_heterozygosity

    # # calculate the mean Fst within the chromosome
    # mean_fst=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_fst)
    # # calculate the median Fst within the chromosome (used '-g' to sort in numerical order while correctly handling scientific notation)
    # median_fst=$(awk '{print $1}' $out_fst | sort -g | awk ' { a[i++]=$1; } END { if (i % 2 == 0) print (a[i/2-1] + a[i/2]) / 2; else print a[int(i/2)]; }')
    # # print a header to the summary file
    # if [ ! -f $summary ]; then
    #     echo -e "chr\tmean_fst\tmedian_fst" > $summary
    # fi
    # echo -e "$chr\t$mean_fst\t$median_fst" >> $summary
done