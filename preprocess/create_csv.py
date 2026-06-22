from os import DirEntry, path, scandir
from pandas import Series, DataFrame, concat, read_csv, to_numeric
from pandas.core.groupby import DataFrameGroupBy 

CURRENT_DIR: str = path.dirname(path.abspath(__file__))
ROOT_DIR: str = path.dirname(CURRENT_DIR)
INPUT_DIR: str = path.join(ROOT_DIR, "data", "output")
OUTPUT_DIR: str = path.join(ROOT_DIR, "data")
FASTA_FILE_PATH: str = path.join(ROOT_DIR, "data", "input", "all_sequences_clean.fasta")

COLUMN_NAMES: list[str] = [
    "Sequence_ID", "Start_X", "End_X", "Start_Y", "End_Y", "Length", "Repeat", "Complement"
]

ALL_PREFIX_FILTER: str = "all_"
COMPLETE_PREFIX_FILTER: str = "complete_"
ALL_SEQUENCES_CSV_FILE_NAME: str = "dataset_all_sequences.csv"
COMPLETE_SEQUENCES_CSV_FILE_NAME: str = "dataset_complete_sequences.csv"

def parse_fasta_species_mapping(fasta_path: str) -> dict[str, str]:
    """Parses the clean FASTA file headers to dynamically build a mapping dictionary
    between Sequence IDs and their official virus names.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        dict[str, str]: A dictionary mapping Sequence_ID to Virus_Name.

    Raises:
        FileNotFoundError: If the specified FASTA file does not exist.    
        Exception: For any other unexpected errors during file reading or parsing.
    """
    mapping_dict: dict[str, str] = {}
    
    if not path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    try:
        with open(fasta_path, "r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                if line.startswith(">"):
                    # Header format: >NC_077666.1 |Puumala virus CG1820 virus M genome segment...
                    header_content: str = line.lstrip(">").strip()
                    if " |" in header_content:
                        parts: list[str] = header_content.split(" |", 1)
                        seq_id: str = parts[0].strip()
                        desc_lower: str = parts[1].strip().lower()
                        
                        if "puumala" in desc_lower:
                            mapping_dict[seq_id] = "Orthohantavirus puumalaense"
                        elif "hantaan" in desc_lower:
                            mapping_dict[seq_id] = "Orthohantavirus hantanense"
                        elif "dobrava" in desc_lower:
                            mapping_dict[seq_id] = "Orthohantavirus dobravaense"
                        elif "sin nombre" in desc_lower or "sinnombre" in desc_lower:
                            mapping_dict[seq_id] = "Orthohantavirus sinnombreense"
                        else:
                            # Safe biological fallback using the first two words if an unrecognized virus appears
                            words: list[str] = desc_lower.split()
                            mapping_dict[seq_id] = " ".join(words[:2]) if len(words) >= 2 else desc_lower
                            
        return mapping_dict
    except Exception as e:
        raise e


def process_single_analysis_file(file_path: str) -> DataFrame:
    """Loads a single raw StatRepeats output file, bypasses headers dynamically,
    and isolates valid comma-separated records.

    Args:
        file_path (str): The path to the text file.

    Returns:
        DataFrame: A cleaned dataframe containing specific repeat rows.
    """
    df: DataFrame = read_csv(
        file_path,
        sep=",",
        header=None,
        names=COLUMN_NAMES,
        on_bad_lines="skip",
        engine="c",
        low_memory=False
    )

    if df.empty:
        return DataFrame()

    numeric_cols: list[str] = ["Start_X", "End_X", "Start_Y", "End_Y", "Length"]
    for col in numeric_cols:
        df[col] = to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Length"])
    df["Length"] = df["Length"].astype(int)
    return df


def extract_features(df: DataFrame, prefix: str) -> DataFrame:
    """Aggregates structural repeat records into unique, sequence-level 
    mathematical features (counts, averages, and length distributions).

    Args:
        df (DataFrame): Cleaned dataframe of individual repeat items.
        prefix (str): Prefix identifying the analysis type (dn, dc, in, ic).

    Returns:
        DataFrame: Aggregated statistical matrix indexed by Sequence_ID.
    """
    grouped: DataFrameGroupBy = df.groupby("Sequence_ID")

    total_count: Series[int] = grouped.size()
    avg_length: Series[float] = grouped["Length"].mean()
    max_length: Series[int] = grouped["Length"].max()

    # Generate the structural feature space: frequency distribution of specific lengths
    pivot_df: DataFrame = df.pivot_table(
        index="Sequence_ID", 
        columns="Length", 
        aggfunc="size", 
        fill_value=0
    )

    if isinstance(pivot_df, Series):
        pivot_df = pivot_df.to_frame()

    pivot_df = pivot_df.add_prefix(f"{prefix}_len_")

    feature_dict: dict[str, Series] = {
        f"{prefix}_total_count": total_count,
        f"{prefix}_avg_length": avg_length,
        f"{prefix}_max_length": max_length
    }

    base_features_df: DataFrame = DataFrame(feature_dict)
    features: DataFrame = concat([base_features_df, pivot_df], axis=1)
    return features


def build_dataset(file_prefix_filter: str, species_map: dict[str, str]) -> DataFrame:
    """Traverses the output directory via fast system scanners to align all 4 cross-sectional
    analyses (dn, dc, in, ic) into a unified classification grid for a specific prefix.

    Args:
        file_prefix_filter (str): Filename prefix filter ('all_' or 'complete_').
        species_map (dict[str, str]): Dynamic map built from FASTA headers.

    Returns:
        DataFrame: Completed analytical dataset.
    
    Raises:
        Exception: If any unexpected error occurs during file processing or data aggregation.
    """
    feature_matrices: list[DataFrame] = []
    analysis_types: list[str] = ["dn", "dc", "in", "ic"]

    try:
        entry: DirEntry
        with scandir(INPUT_DIR) as entries:
            for entry in entries:
                entry_name_lower: str = str(entry.name).lower()

                if (
                    entry.is_file() 
                    and entry_name_lower.startswith(file_prefix_filter) 
                    and entry_name_lower.endswith(".txt")
                ):
                    parts: list[str] = entry_name_lower.split("_")
                    if len(parts) < 2:
                        continue

                    analysis_type: str = parts[1].replace(".txt", "")
                    if analysis_type not in analysis_types:
                        continue
                    
                    df: DataFrame = process_single_analysis_file(entry.path)
                    if df.empty:
                        print(f"WARNING: DataFrame is empty for file: {entry_name_lower}")
                        continue

                    aggregated_features: DataFrame = extract_features(df, analysis_type)
                    feature_matrices.append(aggregated_features)

        if not feature_matrices:
            print(f"WARNING: No feature matrices found for prefix: {file_prefix_filter}")
            return DataFrame()

        # Combines all feature matrices into a single DataFrame, ensuring alignment on Sequence_ID
        final_df: DataFrame = feature_matrices[0]
        for next_matrix in feature_matrices[1:]:
            final_df = final_df.join(next_matrix, how="outer")

        final_df = final_df.fillna(0)

        final_df = final_df.reset_index()
        final_df["Virus_Species"] = final_df["Sequence_ID"].map(species_map).fillna("Unknown_Orthohantavirus")
        return final_df
        
    except Exception as e:
        raise e


def main() -> None:
    """Main execution function to generate two comprehensive CSV datasets from StatRepeats outputs.

    Returns:
        None
    """
    if not path.exists(INPUT_DIR):
        print(f"Error: Input directory {INPUT_DIR} does not exist.")
        return

    print("Parsing dynamic species mapping from FASTA...")
    species_map: dict[str, str] = parse_fasta_species_mapping(FASTA_FILE_PATH)
    print(f"Successfully mapped {len(species_map)} unique sequences from FASTA file headers.")

    print("Generating Dataset 1: All Sequences...")
    all_sequences_dataset: DataFrame = build_dataset(file_prefix_filter=ALL_PREFIX_FILTER, species_map=species_map)
    if not all_sequences_dataset.empty:
        all_sequences_dataset.to_csv(path.join(OUTPUT_DIR, ALL_SEQUENCES_CSV_FILE_NAME), index=False)
        print("Dataset 1 successfully saved.\n")
    else:
        print("Warning: Dataset 1 is empty. No CSV file was generated.\n")

    print("Generating Dataset 2: Complete Sequences Only...")
    complete_sequences_dataset: DataFrame = build_dataset(file_prefix_filter=COMPLETE_PREFIX_FILTER, species_map=species_map)
    if not complete_sequences_dataset.empty:
        complete_sequences_dataset.to_csv(path.join(OUTPUT_DIR, COMPLETE_SEQUENCES_CSV_FILE_NAME), index=False)
        print("Dataset 2 successfully saved.")
    else:
        print("Warning: Dataset 2 is empty. No CSV file was generated.")


if __name__ == "__main__":
    main()