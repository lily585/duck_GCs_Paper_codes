#!/bin/bash

# Define output directories
output_dir="./tissue_specific/specific_state"
merged_dir="./tissue_specific/merged_results"

# Create output directories
mkdir -p "$output_dir"
mkdir -p "$merged_dir"

# Read organization list from file
mapfile -t organizations < list

# Process each organization in the list
for org in "${organizations[@]}"; do
    echo "Processing organization: $org"
    
    # Process each chromatin state (1-15)
    for state in {1..15}; do
        # Define current organization's state file
        current_file="./tissue_specific/tissue_split_state/${org}_15_segments.bed_F${state}.txt"
        
        # Define merged state file for other organizations
        merged_file="$merged_dir/merged_state_${state}.txt"
        > "$merged_file"  # Clear previous content
        
        # Create merged state file from other organizations
        for other_org in "${organizations[@]}"; do
            # Skip current organization
            if [[ "$other_org" != "$org" ]]; then
                other_file="./tissue_specific/tissue_split_state/${other_org}_15_segments.bed_F${state}.txt"
                if [[ -f "$other_file" ]]; then
                    # Append other organization's state data
                    cat "$other_file" >> "$merged_file"
                fi
            fi
        done
        
        # Sort and deduplicate merged state file
        sort -u "$merged_file" -o "$merged_file"
        
        # Skip if current file doesn't exist
        if [[ ! -f "$current_file" ]]; then
            echo "File $current_file does not exist, skipping."
            continue
        fi
        
        # Find tissue-specific regions (non-overlapping with others)
        output_file="$output_dir/non_overlapping_${org}_state_${state}.txt"
        bedtools intersect -a "$current_file" -b "$merged_file" -v > "$output_file"
        
        echo "Created tissue-specific regions: $output_file"
    done
done

echo "All processing completed. Results in $output_dir, merged files in $merged_dir."
