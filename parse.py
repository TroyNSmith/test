import os
import re

# Constants
EH_TO_KCAL = 627.509474

def parse_file(filepath, pattern):
    """Helper to find a pattern and return the last match."""
    results = []
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        for line in f:
            if pattern in line:
                results.append(line.strip())
    return results[-1] if results else None

def process_folders():
    base_dir = 'results/'
    
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} not found.")
        return

    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        
        # Only process directories that don't have result.txt
        if os.path.isdir(folder_path) and not os.path.exists(os.path.join(folder_path, 'result.txt')):
            print(f"Processing: {folder}")
            
            # 1. Parse SP Energy (calc.log)
            sp_line = parse_file(os.path.join(folder_path, 'calc.log'), "FINAL SINGLE POINT ENERGY")
            # 2. Parse ZPE (freq.log)
            zpe_line = parse_file(os.path.join(folder_path, 'freq.log'), "Zero point energy")
            
            if sp_line and zpe_line:
                # Extract values using regex
                sp_eh = float(re.findall(r"[-+]?\d*\.\d+|\d+", sp_line)[0])
                # We take the Eh value (first number) from the ZPE line
                zpe_eh = float(re.findall(r"[-+]?\d*\.\d+|\d+", zpe_line)[0])
                
                sp_kcal = sp_eh * EH_TO_KCAL
                zpe_kcal = zpe_eh * EH_TO_KCAL
                total_energy = sp_kcal + zpe_kcal
                
                output_lines = [
                    f"SP Energy: {sp_eh} Eh ({sp_kcal:.4f} kcal/mol)",
                    f"ZPE: {zpe_eh} Eh ({zpe_kcal:.4f} kcal/mol)",
                    f"Total Energy: {total_energy:.4f} kcal/mol"
                ]
                
                # 3. Handle Transition States (Imaginary Modes)
                if folder.startswith('T'):
                    freq_path = os.path.join(folder_path, 'freq.log')
                    with open(freq_path, 'r') as f:
                        for line in f:
                            if "***imaginary mode***" in line:
                                output_lines.append(f"Imaginary Frequency: {line.strip()}")
                
                # Write to result.txt
                with open(os.path.join(folder_path, 'result.txt'), 'w') as f_out:
                    f_out.write("\n".join(output_lines))
            else:
                print(f"Skipping {folder}: Missing required log data.")

if __name__ == "__main__":
    process_folders()