import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

query_fasta_path = sys.argv[1] # e.g. "AF-Q9Y6V0-F6-model_v1.pdb.gz.pdb.to_align.fa"
blast_file_path = sys.argv[2] # e.g. "blast-out.b6"
homo_fasta_path = sys.argv[3] # e.g. "Homo.fasta"
output_path_for_msa = sys.argv[4] # e.g. "AF-Q9Y6V0-F6-model_v1.pdb.gz.pdb.to_msa.fa"

df = pd.read_table(blast_file_path,header=None)
cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'.strip().split(' ')
df.columns = cols

df["qseqid2"] = df["qseqid"].str.split('-').str[1]

df2 = df.iloc[1:]
df4 = df2[df2['evalue']<1E-62] #Minimum E-value

df5 = df4[np.invert(df4["qseqid2"]==df4["sseqid"])] #Get rid of identical matches to query
df6 = pd.DataFrame(df5["sseqid"].drop_duplicates())

seqList = list()
for record in SeqIO.parse(homo_fasta_path, "fasta"):
    seqList.append((record.id,record.description,record))

df_db = pd.DataFrame(seqList)
df_db["sseqid"] = df_db[0]

df7 = df6.merge(df_db,on="sseqid")

query_fasta = SeqIO.parse(query_fasta_path, "fasta")
listy = list()
for record in query_fasta:
    listy.append(record)

query_record = listy[0]

write_me_list = df7[2].tolist()

record_list = [query_record,] + write_me_list

with open(output_path_for_msa, "w") as output_handle:
    SeqIO.write(record_list, output_handle, "fasta")
