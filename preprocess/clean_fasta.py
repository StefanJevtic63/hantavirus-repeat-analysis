from os import path
from typing import List, Set

INPUT_DIR: str = "data/input"


def clean_fasta_to_strict_alphabet(input_path: str, output_path: str) -> None:
    """Reads a FASTA file and transforms it into a strict single-line sequence format.
    Enforces a valid nucleotide alphabet (A, C, G, T, U) and eliminates structural 
    anomalies that trigger runtime errors and buffering issues in the C++ parser.

    Args:
        input_path (str): The path to the source FASTA file.
        output_path (str): The path where the strictly cleaned FASTA file will be saved.

    Returns:
        None

    Raises:
        FileNotFoundError: If the input_path does not exist.
    """
    if not path.exists(input_path):
        raise FileNotFoundError(f"The source file was not found: {input_path}")

    cleaned_lines: List[str] = []
    current_sequence_chunks: List[str] = []
    valid_bases: Set[str] = {"A", "C", "G", "T", "U"}
    
    with open(input_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            stripped_line: str = line.strip()
            
            if not stripped_line:
                continue
                
            if stripped_line.startswith(">"):
                # If a previous sequence exists, merge its chunks and finalize the entry
                if current_sequence_chunks:
                    cleaned_lines.append("".join(current_sequence_chunks) + "\n")
                    current_sequence_chunks.clear()
                cleaned_lines.append(stripped_line + "\n")
            else:
                upper_line: str = stripped_line.upper()
                sanitized_chunk: str = "".join(char for char in upper_line if char in valid_bases)
                if sanitized_chunk:
                    current_sequence_chunks.append(sanitized_chunk)

        # Append the final sequence entry remaining in the buffer
        if current_sequence_chunks:
            cleaned_lines.append("".join(current_sequence_chunks))

    # Eliminate any trailing newline character at the end of the file buffer
    if cleaned_lines and cleaned_lines[-1].endswith("\n"):
        cleaned_lines[-1] = cleaned_lines[-1].rstrip("\r\n")

    with open(output_path, "w", encoding="utf-8", newline="\n") as file:
        file.writelines(cleaned_lines)
    
    print(f"Cleaned FASTA file successfully saved to: {output_path}")


def main() -> None:
    """Main execution function to perform strict alphabet preprocessing on hantavirus datasets.
    
    Returns:
        None
    """    
    all_sequences_input: str = path.join(INPUT_DIR, "all_sequences.fasta")
    all_sequences_output: str = path.join(INPUT_DIR, "all_sequences_clean.fasta")
    
    complete_sequences_input: str = path.join(INPUT_DIR, "complete_sequences.fasta")
    complete_sequences_output: str = path.join(INPUT_DIR, "complete_sequences_clean.fasta")

    try:
        clean_fasta_to_strict_alphabet(all_sequences_input, all_sequences_output)
        clean_fasta_to_strict_alphabet(complete_sequences_input, complete_sequences_output)
    except Exception as error:
        print(f"An error occurred during strict preprocessing: {error}")


if __name__ == "__main__":
    main()