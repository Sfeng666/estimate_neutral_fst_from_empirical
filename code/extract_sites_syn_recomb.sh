#!/bin/bash

# Define input and output directories
input_dir="../data/input"
output_dir="../data/output"
count_table_prefix="AllCounts_FR_NoInv_"
reference_bed="$input_dir/all_sites_chr_dm6.bed"
syn_bed="$output_dir/all_sites_chr_dm6_rmlowrecomb_syn.bed"
sample_size=18  # Define your sample size

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Temporary file to store combined output
combined_output="$output_dir/AllCounts_FR_NoInv_combined_rmlowrecomb_syn.txt"
> $combined_output

# Read positions from syn_bed
declare -A pos_rmrecomb_syn
while read -r chr pos; do
    pos_rmrecomb_syn["$chr"]+="$pos "
done < <(awk '{print $1, $2}' $syn_bed)

# Read positions from reference_bed
declare -A pos_ct_tab
while read -r chr start end; do
    for ((i=start; i<end; i++)); do
        pos_ct_tab["$chr"]+="$i "
    done
done < <(awk '{print $1, $2, $3}' $reference_bed)

# Function to extract 4-fold synonymous sites for a given chromosome
extract_sites() {
    local chr=$1
    local count_table="$input_dir/${count_table_prefix}${chr}.txt"
    local output_table="$output_dir/${count_table_prefix}${chr}_rmlowrecomb_syn.txt"

    # Extract rows of count table whose genomic positions match those of 4-fold synonymous sites with low recombination rates
    local i=0
    while read -r line; do
        if [[ -n "$line" ]]; then
            local genomic_position=$(echo ${pos_ct_tab["$chr"]} | cut -d' ' -f$((i+1)))
            i=$((i+1))
            if [[ " ${pos_rmrecomb_syn["$chr"]} " =~ " $genomic_position " ]]; then
                # Also filter for sites that have a total allele count >= sample_size
                local allele_counts=($(echo $line | tr ' ' '\n'))
                local total_count=$(IFS=+; echo "$((${allele_counts[*]}))")
                if [[ $total_count -eq $sample_size ]]; then
                    echo "$line" >> "$output_table"
                elif [[ $total_count -gt $sample_size ]]; then
                    local alleles=()
                    for ((j=0; j<${#allele_counts[@]}; j++)); do
                        for ((k=0; k<${allele_counts[$j]}; k++)); do
                            alleles+=($j)
                        done
                    done
                    local downsampled_alleles=($(shuf -n $sample_size -e "${alleles[@]}"))
                    local downsampled_counts=($(for ((j=0; j<${#allele_counts[@]}; j++)); do echo 0; done))
                    for allele in "${downsampled_alleles[@]}"; do
                        downsampled_counts[$allele]=$((downsampled_counts[$allele]+1))
                    done
                    echo "${downsampled_counts[*]}" | tr ' ' '\t' >> "$output_table"
                fi
            fi
        fi
    done < "$count_table"

    # Append the per-chromosome output to a combined output
    cat "$output_table" >> "$combined_output"
}

# Extract sites for each chromosome
chromosomes=$(awk '{print $1}' $reference_bed | sort | uniq)
for chr in $chromosomes; do
    extract_sites $chr
done

echo "Extraction complete. Output files are located in $output_dir."