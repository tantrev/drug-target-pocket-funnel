import pandas as pd
from Bio import PDB
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress BioPython warnings about PDB format
warnings.simplefilter('ignore', PDBConstructionWarning)

class BFactorModifier:
    def __init__(self, input_csv, input_pdb, output_pdb, residue_offset=0, skipfooter=0, target_chain='A', delimiter=','):
        self.input_csv = input_csv
        self.input_pdb = input_pdb
        self.output_pdb = output_pdb
        self.residue_offset = residue_offset
        self.skipfooter = skipfooter
        self.target_chain = target_chain
        self.delimiter = delimiter
        self.data = None
        self.structure = None
        
    def read_data(self):
        """Read and process data from CSV file."""
        try:
            df = pd.read_csv(
                self.input_csv,
                sep=self.delimiter,  # Use specified delimiter
                names=['residue_number', 'amino_acid_letter', 'value'],
                header=None,
                skipfooter=self.skipfooter,
                engine='python'
            )
            
            # Print first few rows for debugging
            print(f"First few rows of {self.input_csv}:")
            print(df.head())
            
            # Convert residue_number to numeric, forcing errors to NaN
            df['residue_number'] = pd.to_numeric(df['residue_number'], errors='coerce')
            
            # Drop rows where residue_number couldn't be converted
            df = df.dropna(subset=['residue_number'])
            
            # Convert residue_number to integers and apply offset
            df['residue_number'] = df['residue_number'].astype(int) + self.residue_offset
            
            # Create dictionary mapping residue numbers to values
            self.data = dict(zip(df['residue_number'], df['value']))
            print(f"Loaded data for {len(self.data)} residues")
            
        except Exception as e:
            raise RuntimeError(f"Error processing data: {str(e)}")
    
    def read_pdb(self):
        """Read PDB structure using BioPython."""
        try:
            parser = PDB.PDBParser(QUIET=True)
            self.structure = parser.get_structure('protein', self.input_pdb)
        except Exception as e:
            raise RuntimeError(f"Error reading PDB structure: {str(e)}")

    def write_pdb(self):
        """Write modified PDB with explicit B-factor column formatting."""
        try:
            modified_lines = []
            with open(self.input_pdb, 'r') as f:
                for line in f:
                    if line.startswith('ATOM  ') and line[21] == self.target_chain:
                        res_num = int(line[22:26])
                        value = self.data.get(res_num, 0.00)
                        
                        new_line = (
                            line[:60] +  # Keep everything before B-factor
                            f"{value:6.2f}" +  # Format B-factor (6 chars, 2 decimal places)
                            line[66:]  # Keep everything after B-factor
                        )
                        modified_lines.append(new_line)
                    elif not line.startswith(('ANISOU', 'TER', 'END')):
                        modified_lines.append(line)
            
            with open(self.output_pdb, 'w') as f:
                f.writelines(modified_lines)
                f.write('END\n')  # Ensure proper PDB format termination
                
        except Exception as e:
            raise RuntimeError(f"Error writing modified PDB: {str(e)}")
    
    def process(self):
        """Execute the modification pipeline."""
        try:
            print("Reading data...")
            self.read_data()
            
            print("Reading PDB structure...")
            self.read_pdb()
            
            print("Writing modified structure...")
            self.write_pdb()
            
            print(f"Successfully modified PDB file. Output written to: {self.output_pdb}")
            
        except Exception as e:
            print(f"Error during processing: {str(e)}")
            raise

# Processing with different configurations
if __name__ == "__main__":
    
    modifier1 = BFactorModifier(
        input_csv="../outputs/hhpred_querytemplate_il12b_domain_only_filtered_al2co_output.csv.txt",
        input_pdb="../pdb_processed/AF-P29460-F1-model_v1.pdb",
        output_pdb="../outputs/conservation_scored_AF-P29460-F1-model_v1.pdb",
        residue_offset=120,
        skipfooter=16,
        target_chain='A',
        delimiter=r'\s+'
    )
    modifier1.process()