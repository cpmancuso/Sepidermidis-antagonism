import pandas as pd
import os
import csv

def extract_knownclusterblast_hits(base_dir, output_csv):
    """
    Parse .txt files in `knownclusterblast` folders and extract cluster details.
    Save results to a CSV file, including the lowest %identity from BLAST hits.
    """
    results = []

    # Traverse the base directory
    for root, dirs, files in os.walk(base_dir):
        if "knownclusterblast" in root:  # Focus only on `knownclusterblast` folders
            # Extract the sample name (parent directory of `knownclusterblast`)
            sample_name = os.path.basename(os.path.dirname(root))

            for file in files:
                if file.endswith(".txt"):  # Process only .txt files
                    file_path = os.path.join(root, file)

                    # Open and parse the file for significant hits
                    with open(file_path, 'r') as f:
                        processed_regions = set()  # Set to track processed contigs
                        lines = f.readlines()
                        parsing = False
                        cluster_details = {}
                        lowest_identity = None  # Variable to store the lowest %identity
                        
                        # Flag to determine if we are in the BLAST hits table
                        in_blast_hits_table = False

                        for line in lines:


                            # Process the rest of the file for cluster details
                            if not parsing and "Significant hits:" in line:
                                next_line_index = lines.index(line) + 1
                                if next_line_index < len(lines) and lines[next_line_index].strip():
                                    parsing = True
                                else:
                                    break

                            if parsing and ("Cumulative BLAST Score" not in cluster_details):
                                # Extract details from subsequent lines
                                if line.strip().startswith("Source:"):
                                    cluster_details["Source"] = line.split(":", 1)[1].strip()
                                elif line.strip().startswith("Type:"):
                                    cluster_details["Type"] = line.split(":", 1)[1].strip()
                                elif line.strip().startswith("Number of proteins with BLAST hits to this cluster:"):
                                    cluster_details["Proteins with BLAST Hits"] = line.split(":", 1)[1].strip()
                                elif line.strip().startswith("Cumulative BLAST score:"):
                                    cluster_details["Cumulative BLAST Score"] = line.split(":", 1)[1].strip()

                            # Check if we have reached the BLAST hits table
                            if line.strip().startswith("Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):"):
                                in_blast_hits_table = True
                                continue  # Skip this header line

                            if in_blast_hits_table and ("Lowest Identity" not in cluster_details):
                                # Stop parsing the BLAST hits table once we hit a blank line or another section
                                if not line.strip():
                                    in_blast_hits_table = False
                                    # Add the lowest %identity to the cluster details if it was found
                                    if (lowest_identity is not None) and ("Lowest Identity" not in cluster_details):
                                        cluster_details["Lowest Identity"] = lowest_identity

                                # Parse the line of BLAST hits (tab-delimited)
                                columns = line.split("\t")
                                if len(columns) >= 6:
                                    try:
                                        # Extract the %identity (third column) from the BLAST hits table
                                        identity = float(columns[2])
                                        # Update the lowest_identity if necessary
                                        if lowest_identity is None or identity < lowest_identity:
                                            lowest_identity = identity
                                    except ValueError:
                                        pass  # In case of invalid data, continue to the next line

                        # Commit the cluster details if populated
                        if cluster_details:
                            # Extract the contig name from the first line
                            temp = lines[0].strip()
                            contig_name = temp[24:]

                            # Extract the region name from filename
                            region_name = os.path.splitext(file)[0]
                            cluster_details["Region Name"] = region_name

                            # Only append the first hit for each contig
                            if (region_name not in processed_regions) and ("Cumulative BLAST Score" in cluster_details):
                                cluster_details["Sample Name"] = sample_name
                                cluster_details["Contig Name"] = contig_name
                                # Extract filename without .txt
                                region_name = os.path.splitext(file)[0]
                                cluster_details["Region Name"] = region_name
                                results.append(cluster_details)
                                processed_regions.add(region_name)  # Mark this region as processed

    # Save results to CSV
    if results:
        fieldnames = [
            "Sample Name",
            "Contig Name",
            "Region Name",
            "Source",
            "Type",
            "Proteins with BLAST Hits",
            "Cumulative BLAST Score",
            "Lowest Identity"
        ]
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"Results saved to {output_csv}")
    else:
        print("No significant hits found.")




def collate_bigscape_results(folder_path, output_tsv):
    """
    Parse .txt files in `bigscape` folders and extract family details.
    Save results to a TSV file.
    """
    base_path = "03-bigscape/output_files/"

    # Get the single subfolder in the base path
    subfolders = [entry for entry in os.listdir(folder_path) if os.path.isdir(os.path.join(base_path, entry))]

    # Check if a subfolder exists
    if subfolders:
        bigscape_folder_path = os.path.join(base_path, subfolders[0])  # Since there's only one, take the first
        print(f"Found folder: {bigscape_folder_path}")
    else:
        print("No subfolder found.")


    # Initialize a list to store data from clustering files
    clustering_data = []

    # Initialize a variable to store the record_annotations table
    record_annotations = None

    # Traverse the directory and its subdirectories
    for root, _, files in os.walk(bigscape_folder_path):
        for file in files:
            file_path = os.path.join(root, file)

            # Identify the record_annotations.tsv file
            if file == "record_annotations.tsv":
                record_annotations = pd.read_csv(file_path, sep='\t')
                print(f"Loaded {file} from {file_path}")

            # Identify clustering files named XX_clustering_c0.3
            elif file.endswith("_clustering_c0.3.tsv"):
                clustering_df = pd.read_csv(file_path, sep='\t')
                clustering_data.append(clustering_df)
                print(f"Added {file} from {file_path}")

    # Combine all clustering data into a single DataFrame
    if clustering_data:
        combined_clustering = pd.concat(clustering_data, ignore_index=True)
        print("Combined all clustering files into one table.")
    else:
        print("No clustering files found.")
        return

    # Merge with record_annotations table if it exists
    if record_annotations is not None:
        # Ensure there's a common key to merge on (update 'common_key' to your actual key column name)
        merged_table = pd.merge(record_annotations, combined_clustering, on="GBK", how="inner")
        print("Merged the combined clustering table with the record_annotations table.")
    else:
        print("record_annotations.tsv not found. Merging skipped.")
        return

    # Save the merged table to a new TSV file

    merged_table.to_csv(output_tsv, sep='\t', index=False)
    print(f"Merged table saved to {output_tsv}")


def collate_by_column_with_lineage(df, group_by_column):
    grouped = df[df[group_by_column].notna()].groupby(group_by_column)
    collated = grouped.agg({
        'Sample_Name': lambda x: ', '.join(x),  # List all sample names
        **{col: 'first' for col in df.columns if col != 'Sample_Name' and col != group_by_column}
    }).reset_index()
    
    # Add Lineage column by extracting and sorting numbers from Sample_Name
    collated['Lineage'] = collated['Sample_Name'].apply(
        lambda x: ', '.join(sorted(
            [name.split('_')[-1] if 'clade' in name else name[2:] for name in x.split(', ')],
            key=int
        ))
    )
    return collated


def main():
    """
    Parse bigscape and knownclusterblast results.
    Merge results by family and by 
    Save results to a TSV file.
    """
    output_path = "05-summarize-results"
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # parse BigScape results
    bigscape_folder_path = "03-bigscape/output_files/"
    bigscape_tsv_file = output_path + "/" + "bigscape_collated_results.tsv"
    collate_bigscape_results(bigscape_folder_path,bigscape_tsv_file)

    # parse antiSMASH results
    knowncluster_folder = "01-antismash"
    knowncluster_csv_file = output_path + "/" + "significant_knowncluster_hits.csv"
    extract_knownclusterblast_hits(knowncluster_folder, knowncluster_csv_file)

    # Load the files
    tsv_data = pd.read_csv(bigscape_tsv_file, sep="\t")
    tsv_data.columns = tsv_data.columns.str.replace(' ', '_')
    csv_data = pd.read_csv(knowncluster_csv_file)
    csv_data.columns = csv_data.columns.str.replace(' ', '_')

    # Parse the Region_Name in the CSV file to extract Contig and Region_Code
    csv_data['Contig'] = csv_data['Region_Name'].str.extract(r"^(contig_\d+)_")
    csv_data['Region_Code'] = csv_data['Region_Name'].str.extract(r"_c(\d+)")

    # Map Region_Code (e.g., c1) to region number (e.g., .region001)
    csv_data['Region'] = csv_data['Region_Code'].apply(lambda x: f".region{int(x):03}")

    # Create the GBK column by combining Sample_Name, Contig, and Region
    csv_data['GBK'] = csv_data['Sample_Name'] + "_" + csv_data['Contig'] + csv_data['Region']

    # Merge the tables on GBK
    merged_data = pd.merge(tsv_data, csv_data, on='GBK', how='outer')

    # Drop the specified columns
    columns_to_drop = [
        'Full_Name_x', 'Record_Type_y', 'Organism', 'Taxonomy', 'Description',
        'Record_Type_x', 'Full_Name_y', 'Record_Number', 'Region_Name', 'Contig_Name',
        'Contig', 'Region_Code', 'Region'
    ]
    merged_data = merged_data.drop(columns=columns_to_drop, errors='ignore')

    # Extract the Sample_Name from GBK
    merged_data['Sample_Name'] = merged_data['GBK'].str.extract(r"^(.*)_contig")

    # Save the merged table
    merged_data_file = output_path + "/" + "merged_tables.csv"
    merged_data.to_csv(merged_data_file, index=False)
    print("Saved as 'merged_tables.csv'")

    # Collate by Family
    collated_by_family = collate_by_column_with_lineage(merged_data, 'Family')
    output_family_file = output_path + "/" + "BGCs_collated_by_family.csv"

    # Collate by Source
    collated_by_source = collate_by_column_with_lineage(merged_data, 'Source')
    output_source_file = output_path + "/" + "BGCs_collated_by_source.csv"

    # Save the resulting tables
    collated_by_family.to_csv(output_family_file, index=False)
    collated_by_source.to_csv(output_source_file, index=False)

    print("Collated tables saved as 'BGCs_collated_by_family.csv' and 'BGCs_collated_by_source.csv'")

if __name__ == '__main__':
    main()