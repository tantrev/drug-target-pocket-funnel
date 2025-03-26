from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read the original alignment
input_file = "../outputs/hhpred_querytemplate_il12b_domain_only_filtered.fasta"
output_file = "../outputs/hhpred_querytemplate_il12b_domain_only_filtered.aln"

# Read alignment
alignment = AlignIO.read(input_file, "fasta")

# Find positions where reference sequence (first sequence) has a residue
ref_seq = str(alignment[0].seq)
keep_positions = [i for i, aa in enumerate(ref_seq) if aa != '-']

# Create new alignment with only those positions
new_records = []
for record in alignment:
    # Extract only the positions we want to keep
    new_seq = ''.join(str(record.seq)[i] for i in keep_positions)
    # Create new record with simplified header
    new_record = SeqRecord(Seq(new_seq), 
                          id=record.id.split()[0],  # Take just the first part of the header
                          description="")
    new_records.append(new_record)

# Create new alignment
new_alignment = MultipleSeqAlignment(new_records)

# Write to CLUSTAL format
AlignIO.write(new_alignment, output_file, "clustal")

print("Conversion complete! Check", output_file)
print(f"Original length: {len(alignment[0].seq)}")
print(f"Trimmed length: {len(new_alignment[0].seq)}")