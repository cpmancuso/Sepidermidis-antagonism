#!/bin/bash

# Source directory
source_dir="known_cluster_hits"

# String to match in the last line
string_to_match="Details:"

# Function to delete text files with the specified string in the last line
delete_files_with_last_line() {
    local dir="$1"
    local match="$2"

    # Delete text files with the specified string in the last line
    find "$dir" -type f -name "*.txt" -exec sh -c '
        for file do
            last_line=$(tail -n 1 "$file")
            if [ "$last_line" = "$1" ]; then
                rm "$file"
            fi
        done
    ' sh "$match" {} +
}

# Call the function to delete files with the specified string in the last line
delete_files_with_last_line "$source_dir" "$string_to_match"