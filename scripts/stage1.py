import pandas as pd
import numpy as np

# Load DGIdb interactions and perform initial filtering
dgidb_interactions = pd.read_table("../data/DGIdb_interactions_12_2023.tsv")
dgidb_interactions["Gene name"] = dgidb_interactions["gene_claim_name"]
dgidb_interactions = dgidb_interactions.dropna(subset=["Gene name"])
inhibitor_interactions = dgidb_interactions[dgidb_interactions["interaction_type"] == "inhibitor"]
inhibitor_interactions = inhibitor_interactions[inhibitor_interactions["approved"]==True]
print("DGIdb filter - unique genes remaining:", inhibitor_interactions["Gene name"].unique().shape[0])

# Load UniProt data and merge with inhibitor interactions
uniprot_data = pd.read_table("../data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.08.31-21.45.11.90.tsv")
uniprot_data["Gene name"] = uniprot_data["Gene Names (primary)"]
uniprot_merged = uniprot_data.merge(inhibitor_interactions, on="Gene name") #inhibitor_interactions
print("UniProt structure filter - unique genes remaining:", uniprot_merged["Gene name"].unique().shape[0])

#Antibody structure filter
thera_df = pd.read_csv("../data/TheraSAbDab_SeqStruc_OnlineDownload.csv")
si_columns = [col for col in thera_df.columns if "SI Structure" in col]
sub_thera_df = thera_df.loc[thera_df[si_columns].notna().any(axis=1)]
delimiters = r'[\/;,]'
sub_thera_df['Gene name'] = sub_thera_df['Target'].str.split(delimiters, regex=True)
sub_thera_df['Gene name'] = sub_thera_df['Gene name'].apply(lambda lst: [x.strip() for x in lst if x and isinstance(x, str)])
sub_thera_df_exploded = sub_thera_df.explode('Gene name')
anti_filtered = uniprot_merged.merge(sub_thera_df_exploded, on="Gene name")
print("Antibody structure filter - unique genes remaining:", anti_filtered["Gene name"].unique().shape[0])

# Load Probe Miner data, find targets without any probes, and merge with UniProt data
unique_probes = pd.read_csv("../data/probeminer_2021-06-20_unique_probes.csv")
merged_probes = unique_probes.merge(anti_filtered, on="Entry", indicator=True, how='outer')
right_only_probes = merged_probes[merged_probes["_merge"] == "right_only"]
print("Probe Miner filter - unique genes remaining:", right_only_probes["Gene name"].unique().shape[0])

#Filter out proteins with metals
metals = pd.read_table("../data/uniprot-annotation_(type_metal)-filtered-proteome_UP000005640+AND+orga--.tab")
metals["Gene name"] = metals["Gene names  (primary )"]
merged_no_metals = right_only_probes.drop("_merge", axis=1).merge(metals, on="Gene name", how='outer', indicator=True)
merged_no_metals = merged_no_metals[merged_no_metals["_merge"]=="left_only"]
print("Metals filter - unique genes remaining:", merged_no_metals["Gene name"].unique().shape[0])

#Bioinformatics filter - only looks at proteins where no other human protein has significant primary sequence homology
#Also filter out proteins that are longer than AlphaFold's max length for 1 protein prediction (and thus split the prediction into multiple files)
bio_df = pd.read_csv("../outputs/homology_counts.txt", header=None) #22_with_eval_cut.txt
bio_df["Entry_x"] = bio_df[0].str.split("-").str[1].astype('str')
bio_df["Count"] = bio_df[0].str.split(":").str[-1].astype('int')
bio_df2 = bio_df[bio_df["Count"]==1] #Make sure there are no homologous proteins
bio_df3 =  bio_df2["Entry_x"].value_counts()
bio_df4 =  bio_df3[bio_df3 == 1].reset_index() #Make sure the protein is not longer than AlphaFold's max prediction length
bio_df4.columns = ["Entry_x", "Count2"]
bio_df5 = bio_df4.merge(bio_df2, on="Entry_x")
final_filtered = merged_no_metals.merge(bio_df5, on="Entry_x")
print("Homology filter - unique genes remaining:", final_filtered["Gene name"].unique().shape[0])

# Save final filtered 
final_filtered =  final_filtered.copy()
final_filtered["Entry"] = final_filtered["Entry_x"]
final_filtered[["Gene name", "Entry"]].drop_duplicates().to_csv("../outputs/stage1_filtering_gene_candidates.csv", index=None)