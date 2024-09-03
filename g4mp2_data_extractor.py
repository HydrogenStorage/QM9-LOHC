import os
import re
import pandas as pd

def extract_data_from_log(file_path):
    """
    Extracts the required data from a .log file, unless it contains
    the phrase "imaginary frequencies".

    Args:
    file_path (str): The path to the .log file.

    Returns:
    dict: A dictionary containing the extracted data or None if the file contains imaginary frequencies.
    """
    data = {
        "ID": None,
        "SMILES": None,
        "E(ZPE)": None,             # Added E(ZPE)
        "G4MP2(0 K)": None,
        "G4MP2 Energy": None,
        "G4MP2 Enthalpy": None,
        "G4MP2 Free Energy": None
    }

    # Extract the molecule ID from the file name
    file_name = os.path.basename(file_path)
    match = re.search(r"molecule-(\d+)\.log", file_name)
    if match:
        data["ID"] = match.group(1)

    with open(file_path, 'r') as file:
        file_content = file.read()

        # Debug: Print file name and first 200 characters for checking
        print(f"Processing file: {file_path}")
        print(f"File content (first 200 characters): {file_content[:200]}")

        # Check for the presence of imaginary frequencies
        if "imaginary frequencies" in file_content:
            print(f"Imaginary frequencies found in file: {file_path}")
            return None

        # Extract the values and SMILES string if "G4MP2 Enthalpy" is present
        if "G4MP2 Enthalpy" in file_content:
            # Extract G4MP2 and E(ZPE) values
            g4mp2_pattern = (
                r"E\(ZPE\)=\s+([-\d\.]+)\s+E\(Thermal\)=\s+[-\d\.]+\s+"
                r"E\(QCISD\(T\)\)=\s+[-\d\.]+\s+E\(Empiric\)=\s+[-\d\.]+\s+"
                r"DE\(MP2\)=\s+[-\d\.]+\s+DE\(HF\)=\s+[-\d\.]+\s+"
                r"G4MP2\(0 K\)=\s+([-\d\.]+)\s+G4MP2 Energy=\s+([-\d\.]+)\s+"
                r"G4MP2 Enthalpy=\s+([-\d\.]+)\s+G4MP2 Free Energy=\s+([-\d\.]+)"
            )
            g4mp2_match = re.search(g4mp2_pattern, file_content)
            if g4mp2_match:
                data["E(ZPE)"] = float(g4mp2_match.group(1))
                data["G4MP2(0 K)"] = float(g4mp2_match.group(2))
                data["G4MP2 Energy"] = float(g4mp2_match.group(3))
                data["G4MP2 Enthalpy"] = float(g4mp2_match.group(4))
                data["G4MP2 Free Energy"] = float(g4mp2_match.group(5))
                print(f"Extracted data from file: {file_path}")
            else:
                print(f"G4MP2 pattern not found in file: {file_path}")

            # Extract SMILES string
            smiles_pattern = r"Optimization for\s+(.+)"
            smiles_match = re.search(smiles_pattern, file_content)
            if smiles_match:
                data["SMILES"] = smiles_match.group(1).strip()
            else:
                print(f"SMILES pattern not found in file: {file_path}")

    return data if data["G4MP2 Enthalpy"] is not None else None

def process_directory(directory):
    """
    Processes all .log files in a directory and its subdirectories.

    Args:
    directory (str): The path to the directory.

    Returns:
    list: A list of dictionaries containing the extracted data from each file.
    """
    results = []

    for root, dirs, files in os.walk(directory):
        print(f"Current directory being processed: {root}")
        print(f"Subdirectories in current directory: {dirs}")
        print(f"Files in current directory: {files}")

        for file in files:
            if file.endswith(".log"):
                print(f"Found .log file: {file}")
                file_path = os.path.join(root, file)
                data = extract_data_from_log(file_path)
                if data is not None:
                    results.append(data)
                else:
                    print(f"No data extracted from file: {file_path}")
            else:
                print(f"Skipping non-.log file: {file}")

    return results

# Example usage
directory_path = '.'
extracted_data = process_directory(directory_path)

# Check if any data was extracted
if extracted_data:
    # Create DataFrame and save to CSV
    df = pd.DataFrame(extracted_data)
    csv_file_path = "g4mp2_results_new.csv"
    df.to_csv(csv_file_path, index=False)

    # Display the DataFrame
    print(df.head())
else:
    print("No valid data extracted from any log files.")

