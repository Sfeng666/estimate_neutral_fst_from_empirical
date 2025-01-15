#! /bin/bash
# This script prepare input data for the whole analysis

# 1.  prep genotype count tables and genomic coordinates of chromosomes from lab computer
populations=(ZI FR)
for pop in ${populations[@]}; do
    name_pattern="AllCounts_${pop}_NoInv_Chr*.txt"  # pre-defined pattern to match the count data tables
    list_ct_tab=../data/input/source_path2count_table_${pop}.txt  # list of original paths to count data tables for the population
    find /raid10/backups/genepool/DPGP2plus/ -type f -name $name_pattern > $list_ct_tab 2>/dev/null # find the count data tables and save the paths to the list file
    # deduplicate the list of count data tables by chromosome name
    temp_file=$(mktemp) # Temporary file to store deduplicated list
    ## Read the file and deduplicate by chromosome
    awk '
    {
        # Extract the chromosome information (e.g., Chr3R)
        match($0, /Chr[0-9A-Za-z]+/, arr)
        chr = arr[0]
        
        # Keep only the first occurrence of each chromosome
        if (!seen[chr]++) {
            print $0
        }
    }
    ' "$list_ct_tab" > "$temp_file"
    mv "$temp_file" "$list_ct_tab"  # overwrite the original file

    # copy the count data tables to the input data folder
    # also add each site within the chromosome as an entry in a BED file for later position calculation
    bed_all_sites_chr=../data/input/all_sites_chr_dm3.bed
    if [[ -f $bed_all_sites_chr ]]; then
        create_bed=false   # Remove the existing BED file, since both population will generate the same file
    fi

    while read -r ct_tab_path; do
        ct_tab_name=$(basename $ct_tab_path)
        cp $ct_tab_path ../data/input/$ct_tab_name
        if $create_bed; then # Only generate the BED file for one time, since both population will generate the same file
            chr_name=$(echo $ct_tab_name | grep -o 'Chr[0-9A-Za-z]*' | sed 's/Chr/chr/')  # Extract the chromosome name from the file name (assuming the format includes the chromosome name)
            file_length=$(wc -l < $ct_tab_path)   # Calculate the file length (number of lines)
            # generate a BED file with every site in the chromosome as an entry
            for ((i=0; i<file_length; i++)); do
                echo -e "${chr_name}\t${i}\t$((i+1))\t${i}" >> $bed_all_sites_chr # 4th colume is the site position in the chromosome, starting from 0
            done
        fi
    done < $list_ct_tab
done

# 2. download the chain file for liftingover site positions from Drosophila melanogaster ref genome dm3 to dm6
over_chain_file=../data/input/dm3ToDm6.over.chain.gz
url_to_chain_file=http://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
wget $url_to_chain_file -O $over_chain_file

# 3. Liftover both the BED coordinates of full chromosomes and low-recombination regions
bed_all_sites_chr_converted=../data/input/all_sites_chr_dm6.bed
liftOver $bed_all_sites_chr $over_chain_file $bed_all_sites_chr_converted $bed_all_sites_chr_converted.unMapped &
# remove unMapped file if it is empty
if [[ ! -s $bed_all_sites_chr_converted.unMapped ]]; then
    rm $bed_all_sites_chr_converted.unMapped
fi

bed_low_recomb=../data/input/low_recomb_regions_dm3.bed
bed_low_recomb_converted=../data/input/low_recomb_regions_dm6.bed
liftOver $bed_low_recomb $over_chain_file $bed_low_recomb_converted $bed_low_recomb_converted.unMapped
# remove unMapped file if it is empty
if [[ ! -s $bed_low_recomb_converted.unMapped ]]; then
    rm $bed_low_recomb_converted.unMapped
fi

# 4. prepare the annotation sequence file of Drosophila melanogaster reference genome dm6 (forward strand)
annotation_seq="../data/input/for_reannotate_dm6.fa"
conversion_table="../data/input/chr2acc_dm6"
annotation_seq_converted="../data/input/for_reannotate_dm6_converted.fa"

# copy the original annotation sequence file to input folder
cp /home/siyuan/jobs/Dmel_cold_adapt_rna_edit/genomic_annotation_noncoding/output/for_reannotate.fa $annotation_seq
# since the annotation sequence file used a NCBI refseq ID as the header, we need to replace it with the corresponding chromosome name
## first download the chr2acc file from NCBI FTP server
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc -O $conversion_table
## next, convert the NCBI refseq ID to chromosome name in the annotation sequence file
### Create a sed script to perform the replacements
sed_script=$(mktemp)
grep -v '^#' $conversion_table | while IFS=$'\t' read -r chr acc; do
    echo "s/>$acc/>chr$chr/" >> $sed_script
done

### Apply the sed script to the fasta file, and filter out chromosomes not in the conversion table (contains only primary chromosomes)
{
    in_header=false
    while IFS= read -r line; do
        if [[ $line == \>* ]]; then
            acc=${line#>}
            chr=$(grep -P "\t$acc$" $conversion_table | cut -f1)
            if [[ -n $chr ]]; then
                echo ">chr$chr"
                in_header=true
            else
                in_header=false
            fi
        else
            if $in_header; then
                echo $line
            fi
        fi
    done < $annotation_seq
} > $annotation_seq_converted

### Clean up the temporary sed script
rm $sed_script

# 5. Prepare a BED (.bed) file containing genomic positions of 4-fold synonymous sites.
annotation_pos_converted_syn="${annotation_seq_converted%.*}_syn.bed"
awk '/^>/ {chrom=$0; gsub(/^>/, "", chrom); pos=0; next} {
    for (i=1; i<=length($0); i++) {
        if (substr($0, i, 1) == "S")
            print chrom "\t" pos "\t" pos+1
        pos++
    }
}' $annotation_seq_converted > $annotation_pos_converted_syn

# 6. Exclude sites in low-recombination regions, by excluding positions in the BED file of low-recombination regions from the BED file of all sites using `bedtools subtract`.
bed_all_sites_chr_converted_rmlowrecomb="../data/output/$(basename ${bed_all_sites_chr_converted%.*})_rmlowrecomb.bed"
bedtools subtract -a $bed_all_sites_chr_converted -b $bed_low_recomb_converted > $bed_all_sites_chr_converted_rmlowrecomb

# 7. Subset 4-fold synonymous sites for Fst estimation, by extracting positions in the BED file of 4-fold synonymous sites that overlap with the BED output of last step using `bedtools intersect`.
bed_all_sites_chr_converted_rmlowrecomb_syn="${bed_all_sites_chr_converted_rmlowrecomb%.*}_syn.bed"
bedtools intersect -a $bed_all_sites_chr_converted_rmlowrecomb -b $annotation_pos_converted_syn > $bed_all_sites_chr_converted_rmlowrecomb_syn

# optional - check what proportion of sites are analyzed
length_genome=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' /home/siyuan/jobs/Dmel_cold_adapt_rna_edit/genomic_annotation_noncoding/input/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna | awk '{s+=$1} END {print s}')  # calculate the total length of fly gennome (dm6)
num_total_sites_syn=$(wc -l < $annotation_pos_converted_syn) # count the number of all 4-fold synonymous sites in the genome
num_sites_syn_rmlowrecomb=$(wc -l < $bed_all_sites_chr_converted_rmlowrecomb_syn)  # count the number of 4-fold synonymous sites outside the low-recombination regions
num_sites_syn_lowrecomb=$(bedtools intersect -a $bed_low_recomb_converted -b $annotation_pos_converted_syn | wc -l)  # count the number of 4-fold synonymous sites in the low-recombination regions
echo "Proportion of 4-fold synonymous sites across the genome: $(echo "scale=4; $num_total_sites_syn / $length_genome" | bc)"
echo "Proportion of 4-fold synonymous sites in the low-recombination regions: $(echo "scale=4; $num_sites_syn_lowrecomb / $num_total_sites_syn" | bc)"
echo "Proportion of 4-fold synonymous sites outside the low-recombination regions: $(echo "scale=4; $num_sites_syn_rmlowrecomb / $num_total_sites_syn" | bc)"
echo "Proportion of 4-fold synonymous sites in the low-recombination regions compared to the whole genome: $(echo "scale=4; $num_sites_syn_rmlowrecomb / $length_genome" | bc)"
echo "Number of 4-fold synonymous sites in the low-recombination regions: $num_sites_syn_rmlowrecomb"


