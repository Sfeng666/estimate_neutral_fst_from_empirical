#! /bin/bash
# Estimate minor allele frequency spectrum (AFS) for each population

# fixed parameters
pop1="ZI"
pop2="FR"

# demographic parameters that vary across chromosomal arms or combinations of populations
chrs=("chr2L" "chr2R" "chr3L" "chr3R" "chrX" "combined")

for chr in "${chrs[@]}"; do
    in_ct_pop1="../data/output/AllCounts_${pop1}_NoInv_${chr}_rmlowrecomb_syn_seg.txt"
    in_ct_pop2="../data/output/AllCounts_${pop2}_NoInv_${chr}_rmlowrecomb_syn_seg.txt"
    out_afs_pop1="../data/output/afs_${pop1}_${chr}.txt"
    out_afs_pop2="../data/output/afs_${pop2}_${chr}.txt"
    
    # estimate minor allele frequency spectrum (AFS)
    python calc_snp_afs_from_ct.py --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_afs_pop1 $out_afs_pop1 --out_afs_pop2 $out_afs_pop2
done