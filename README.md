# Estimate SNP FST between two populations using nearly-neutral sites
This analysis estimate SNP (i.e., single-site) FST (between two focal populations, here France and Zambia) for 4-fold synonymous sites in moderate/high-recombination regions. 


It serves as an essential step towards generating realistic (adaptive and neutral) QST values from the estimated genome-wide neutral FST distribution, then perform an ROC (Receiver Operating Characteristic) analysis to evaluate how a given estimation model can distinguish between neutral or adaptive (i.e., true outlier) QST.

## Input & Output
Output
* A [summary file](data/output/estimate_fst.report) contains mean FST across all 4-fold synonymous sites excluding low-recombination regions
* For each chromosomal arm, there is an one-column table of estimated SNP-FST ([example](data/output/fst_chr2L.txt)), where each row is a site

Input:
1. Genotype count tables for all chromosomes (one file per chromosome) for each population. Count tables from two populations (e.g., [France](data/input/path2countdata_FR.txt) and [Zambia](data/input/path2countdata_ZI.txt)) are required.
    * file format: each column contains the allele count of A, T, C, G in a population; each row contains the genotype information of one nucleotide, ordered from position 1 to the length of chromosome.
    * pre-filters: genomic regions of specifc genomes of the following features were already excluded from allele count: inversion, IBD/relatedness, apparent heterozygosity, admixture, 3bp filter around called indel.
2. A [BED (.bed) file](data/input/all_sites_chr_dm6.bed) containing coordinates from 0 to the chromosome length of chromosomes that are involved in genotype count tables.
    * Liftover BED coordinates to the same genome assembly version as the annotation sequence file ([code](code/prep_neutral_sites.sh#L49-L54)), if the [original BED file](data/input/all_sites_chr_dm3.bed) is based on a different assembly version.
3. A [BED (.bed) file](data/input/low_recomb_regions_dm6.bed) containing coordinates for low-recombination regions (e.g., low crossing-over regions in [Pool et al., 2017](https://academic.oup.com/mbe/article/34/2/349/2680805)) that will be excluded from the pool of SNPs for FST calculation.
    * Lift over BED coordinates to the same genome assembly version as the annotation sequence file ([code](code/prep_neutral_sites.sh#L56-L62)), if the [original BED files](data/input/low_recomb_regions_dm3.bed) is based on a different assembly version.
4. An [annotation file](data/input/for_reannotate_dm6_converted.fa) of the reference genome sequence (forward strand) in FASTA (.fasta) format, where each nucleotide is replaced by its [annotation symbol](https://github.com/Sfeng666/Dsuz_popgen_GEA/blob/main/genomic_annotation/README.md#this-program-annotates-each-site-of-a-genome-by-nine-pre-defined-categories).
    * (optional) If the header of the [original annotation file](data/input/for_reannotate_dm6.fa) is NCBI refseq ID, you'll need to convert the refseq ID to the corresponding chromosome names using the [NCBI name chromosome name conversion table of the species](data/input/chr2acc_dm6). The code for conversion is in the [input-prep script](code/prep_neutral_sites.sh#L65-L103).
5. A [BED (.bed) file](data/input/for_reannotate_dm6_converted.bed) containing genomic positions of 4-fold synonymous sites. I automatically generated this file using [this shell code](code/prep_neutral_sites.sh#L105-L112).
6. (Optional) An over.chain file to [convert genomic coordinates](https://genome.ucsc.edu/FAQ/FAQdownloads.html#liftOver), if the reference genome version used for genotype count tables is different from that of the annotation file and/or low-recombination region BED file.

Output:

## Method step-by-step
0. Place the required input files into the input folder ([code](code/prep_neutral_sites.sh#L4-L113)).
1. Extract the genomic positions (based on dm6) of 4-fold synonymous sites in moderate-to-high recombination regions ([code](code/prep_neutral_sites.sh#L115-L131))
    1. Exclude sites in low-recombination regions, by excluding positions in the BED file of low-recombination regions from the BED file of all sites using `bedtools subtract`.
    2. Subset 4-fold synonymous sites for Fst estimation, by extracting positions in the BED output of last step that overlaps with positions in the BED file of genome-wide 4-fold synonymous sites using `bedtools intersect`.
2. Extract rows in genotype count tables that represent 4-fold synonymous sites that are not in low-recombination regions ([code](code/extract_sites_syn_recomb.py)).
3. Calculate per-site FST between the two focal populations for each chromosome, and generate mean per-site FST within each chromosome and across chromosomes.

## Note
1. 76.8% 4-fold synonymous sites (based on dm6) locate in low-recombination regions (simply because low-recombination regions occupies a majority of the genome), which are pre-filtered. Only 21.4% 4-fold synonymous sites locate in moderate-to-high recombination regions and are also available in the dm3 genome assembly. Since only 1.2% of the genome are 4-fold synonymous, <u>0.26%</u> of the whole genome are used for estimating SNP FST ([code](code/prep_neutral_sites.sh#L123-L131)).
2. ~0.08% of sites cannot be liftovered from dm3 to dm6. 
3. 183 sites that were on X chromosomes in assembly dm3 were liftover to a weird variant of X chromosome called 'chrX_DS484562v1_random' in dm6. Since the BED file of annotated 4-fold synonymouns sites only contains major chromosomes (in the [name converison file](data/input/chr2acc_dm6)), these 183 sites are not carried over to sites extraction and FST calculation.
4. In step 2 extracting genotype count table, I only keep sites whose total allele count is no less than the pre-defined sample size in both populations. If total allele count exceeds sample size, I downsample the count for each allele proportionally without replacement, so that the downsampled total alelle count = sample_size. 
5. I did not output genomic coordinates with the FST calculations. The coordinates may be more useful when looking for FST outliers beyond neutral regions.

## Environment setup
To set up the environment for this analyses, you could use conda:
```
conda env create -n qst_fst_sim --file env/qst_fst_sim.yml
```

If you have not installed conda, run the following command:
```
# download miniconda
curl -sL \
  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
  "Miniconda3.sh"
```
```  
# install miniconda
bash Miniconda3.sh
```