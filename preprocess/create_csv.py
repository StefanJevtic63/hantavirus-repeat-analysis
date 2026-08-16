from collections import defaultdict
from os import DirEntry, path, scandir

from pandas import DataFrame

CURRENT_DIR: str = path.dirname(path.abspath(__file__))
ROOT_DIR: str = path.dirname(CURRENT_DIR)
INPUT_DIR: str = path.join(ROOT_DIR, "data", "output")
OUTPUT_DIR: str = path.join(ROOT_DIR, "data")

ALL_PREFIX_FILTER: str = "all_"
COMPLETE_PREFIX_FILTER: str = "complete_"

ALL_SEQUENCES_CSV_FILE_NAME: str = "all_sequences.csv"
COMPLETE_SEQUENCES_CSV_FILE_NAME: str = "complete_sequences.csv"


def process_single_analysis_file(file_path: str) -> DataFrame:
    """Loads a single StatRepeats output file and extracts
    (Sequence_ID, Combined_Repeat) pairs.

    Args:
        file_path (str): Path to the StatRepeats output file.

    Returns:
        DataFrame: A DataFrame containing the extracted pairs with columns
                  'Sequence_ID' and 'Repeat'.
    """
    records: list[dict[str, str]] = []

    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            line: str = line.strip()
            if not line or not line[0].isalnum():
                continue

            parts: list[str] = [part.strip() for part in line.split(",")]
            if len(parts) != 8:
                continue

            sequence_id: str = parts[0]
            combined_repeat: str = parts[6] + parts[7]

            records.append(
                {
                    "Sequence_ID": sequence_id,
                    "Repeat": combined_repeat,
                }
            )

    return DataFrame(records)


def build_dataset(file_prefix_filter: str) -> DataFrame:
    """Builds a binary nucleotide repeats feature matrix.

    Rows:
        nucleotide sequences

    Columns:
        all unique combined repeats (LP + RP)

    Values:
        1 -> repeat exists in sequence
        0 -> repeat does not exist
    """
    print("Building dataset for prefix filter:", file_prefix_filter)
    sequence_repeats: dict[str, set[str]] = defaultdict(set) # sequence -> unique repeats
    unique_repeats: set[str] = set() # all unique repeats

    try:
        entry: DirEntry
        with scandir(INPUT_DIR) as entries:
            for entry in entries:
                entry_name = str(entry.name).lower()
                if (
                    not entry.is_file()
                    or not entry_name.startswith(file_prefix_filter)
                    or not entry_name.endswith(".txt")
                ):
                    continue

                df: DataFrame = process_single_analysis_file(entry.path)
                if df.empty:
                    print(f"WARNING: Empty file: {entry.name}")
                    continue

                for _, row in df.iterrows():
                    sequence_id: str = row["Sequence_ID"]
                    repeat: str = row["Repeat"]

                    sequence_repeats[sequence_id].add(repeat)
                    unique_repeats.add(repeat)

                print("Processed:", entry.name, ". Current unique repeats:", len(unique_repeats))

        if not sequence_repeats:
            print(f"WARNING: No data found for prefix '{file_prefix_filter}'.")
            return DataFrame()

        unique_repeat_list: list[str] = sorted(unique_repeats)
        sequence_list: list[str] = sorted(sequence_repeats.keys())
        print(f"Number of unique repeats in '{file_prefix_filter}':", len(unique_repeat_list))
        print(f"Number of sequences in '{file_prefix_filter}':", len(sequence_list))

        dataset: list[dict[str, int]] = []
        for sequence in sequence_list:
            seq_repeats: set[str] = sequence_repeats[sequence]
            dataset.append(
                {
                    repeat: int(repeat in seq_repeats)
                    for repeat in unique_repeat_list
                }
            )

        return DataFrame(dataset, columns=unique_repeat_list)
    except Exception as e:
        raise e

def generate_dataset(csv_file_name: str, prefix_filter: str) -> None:
    """Generates a dataset CSV file based on the given prefix filter.
    
    Args:
        csv_file_name (str): The name of the output CSV file.
        prefix_filter (str): The prefix filter to select input files.
    
    Returns:
        None
    """
    print(f"Generating {csv_file_name}...")
    dataset: DataFrame = build_dataset(prefix_filter)
    if not dataset.empty:
        output_path: str = path.join(OUTPUT_DIR, csv_file_name)
        dataset.to_csv(output_path, index=False)
        print(f"Saved: {output_path}")
    else:
        print(f"No data generated for {csv_file_name}.")

def main() -> None:
    """Generates:
    - all_sequences.csv
    - complete_sequences.csv
    """
    if not path.exists(INPUT_DIR):
        print(f"Error: Input directory '{INPUT_DIR}' does not exist.")
        return

    generate_dataset(ALL_SEQUENCES_CSV_FILE_NAME, ALL_PREFIX_FILTER)
    generate_dataset(COMPLETE_SEQUENCES_CSV_FILE_NAME, COMPLETE_PREFIX_FILTER)


if __name__ == "__main__":
    main()