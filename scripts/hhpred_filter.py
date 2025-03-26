import re
import pandas as pd

def normalize_name(name):
    # Remove isoform indicators but keep identifying numbers
    name = re.sub(r'\sisoform\s+(?:X?\d+|[A-Za-z]+)', '', name)
    name = re.sub(r'\stype\s+\d+', '', name)
    name = re.sub(r'\s+\([^)]+\)', '', name)
    name = re.sub(r'\s+precursor$', '', name)
    name = re.sub(r'\s+X\d+$', '', name)
    return name.strip()

def clean_header(header):
    name_match = re.search(r'>\S+\s+(?:PREDICTED:\s+)?(.*?)(?:\s+\[|$)', header)
    if name_match:
        base_name = name_match.group(1).strip()
        return normalize_name(base_name)
    return header.strip('>')

def process_fasta(fasta_file):
    entries = []
    current_header = ''
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_header:
                    entries.append({
                        'full_header': current_header,
                        'sequence': ''.join(current_seq)
                    })
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
    
    if current_header:
        entries.append({
            'full_header': current_header,
            'sequence': ''.join(current_seq)
        })
    
    df = pd.DataFrame(entries)
    df['clean_name'] = df['full_header'].apply(clean_header)
    filtered_df = df.drop_duplicates(subset=['clean_name'], keep='first')
    
    return df, filtered_df

def save_filtered_fasta(filtered_df, output_file):
    with open(output_file, 'w') as f:
        for _, row in filtered_df.iterrows():
            f.write(f"{row['full_header']}\n")
            f.write(f"{row['sequence']}\n\n")  # Add extra newline between entries

# Process the files
df, filtered_df = process_fasta('../outputs/hhpred_querytemplate_il12b_domain_only.fasta')

# Save filtered results
save_filtered_fasta(filtered_df, '../outputs/hhpred_querytemplate_il12b_domain_only_filtered.fasta')

# Print diagnostic information
print(f"Original entries: {len(df)}")
print(f"Filtered entries: {len(filtered_df)}")
print(f"Unique clean names: {len(filtered_df['clean_name'].unique())}")

# Check for any duplicates that might have slipped through
duplicates = filtered_df[filtered_df['clean_name'].duplicated(keep=False)]
if not duplicates.empty:
    print("\nFound duplicate clean names:")
    for name in duplicates['clean_name'].unique():
        print(f"- {name}")

print("\nSample of clean names:")
print(filtered_df[['clean_name']].head().to_string())

# Additional diagnostic info - check sequence counts
#print("\nChecking input file sequence count:")
#with open('hhpred_querytemplate_il12b_domain_only.fasta', 'r') as f:
    #input_headers = sum(1 for line in f if line.startswith('>'))
#print(f"Number of sequences in input file: {input_headers}")