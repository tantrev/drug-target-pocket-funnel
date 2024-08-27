import pandas as pd
import numpy as np

# Load DGIdb interactions and perform initial filtering
dgidb_interactions = pd.read_table("data/DGIdb_interactions_12_2023.tsv")
dgidb_interactions["Gene name"] = dgidb_interactions["gene_claim_name"]
dgidb_interactions = dgidb_interactions.dropna(subset=["Gene name"])
inhibitor_interactions = dgidb_interactions[dgidb_interactions["interaction_type"] == "inhibitor"]
inhibitor_interactions = inhibitor_interactions[inhibitor_interactions["approved"]==True]
print("DGIdb filter - unique genes remaining:", inhibitor_interactions["Gene name"].unique().shape[0])

# Load UniProt data and merge with inhibitor interactions
uniprot_data = pd.read_table("data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.08.31-21.45.11.90.tsv")
uniprot_data["Gene name"] = uniprot_data["Gene Names (primary)"]
uniprot_merged = uniprot_data.merge(inhibitor_interactions, on="Gene name") #inhibitor_interactions
print("UniProt structure filter - unique genes remaining:", uniprot_merged["Gene name"].unique().shape[0])

# Load Probe Miner data, find targets without any probes, and merge with UniProt data
#probe_miner_data = pd.read_table("data/probeminer_datadump_2021-06-20.txt")
#probe_miner_data["Entry"] = probe_miner_data["UNIPROT_ACCESSION"]
#unique_probes = pd.DataFrame(probe_miner_data["Entry"].drop_duplicates())
#unique_probes.to_csv("probeminer_2021-06-20_unique_probes.txt")
unique_probes = pd.read_csv("data/probeminer_2021-06-20_unique_probes.csv")
merged_probes = unique_probes.merge(uniprot_merged, on="Entry", indicator=True, how='outer')
right_only_probes = merged_probes[merged_probes["_merge"] == "right_only"]
print("Probe Miner filter - unique genes remaining:", right_only_probes["Gene name"].unique().shape[0])

# Load gnomAD data and merge with probes data
gnomad_data = pd.read_table("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
gnomad_data["Gene name"] = gnomad_data["gene"]
gnomad_filtered = gnomad_data[gnomad_data['pLI'] <= .50]
gnomad_merged = right_only_probes.merge(gnomad_filtered, on="Gene name")
print("gnomAD filter - unique genes remaining:", gnomad_merged["Gene name"].unique().shape[0])

# Load and filter OpenTargets data, then merge with gnomAD data
opentargets_data = pd.read_csv("data/BigQuery OpenTarget locus2gene results -12_22_2023.csv")
opentargets_filtered = opentargets_data[opentargets_data["max_y_proba_full_model"] >= 0.5]
opentargets_merged = opentargets_filtered.merge(gnomad_merged, on="gene_id")
print("OpenTargets filter - unique genes remaining:", opentargets_merged["gene_id"].unique().shape[0])

#Filter out proteins with metals
metals = pd.read_table("data/uniprot-annotation_(type_metal)-filtered-proteome_UP000005640+AND+orga--.tab")
metals["Gene name"] = metals["Gene names  (primary )"]
opentargets_merged_no_metals = opentargets_merged.drop("_merge", axis=1).merge(metals, on="Gene name", how='outer', indicator=True)
opentargets_merged_no_metals = opentargets_merged_no_metals[opentargets_merged_no_metals["_merge"]=="left_only"]
print("Metals filter - unique genes remaining:", opentargets_merged_no_metals["Gene name"].unique().shape[0])

# Load mouse map and embryonic lethality data, then perform final merge
mouse_gene_map = pd.read_table("data/HOM_MouseHumanSequence.rpt")
mouse_gene_map["Gene name"] = mouse_gene_map["Symbol"]
embryonic_lethality = pd.read_table("data/embryonic_lethality_MGIhdpQuery_markers_20230416_020549.txt")
embryonic_lethality["Gene name"] = embryonic_lethality["Gene Symbol"]
lethality_merged = embryonic_lethality.merge(mouse_gene_map, on="Gene name")
lethality_final = lethality_merged.merge(mouse_gene_map, on="DB Class Key")
lethality_final["Gene name"] = lethality_final["Gene name_y"]
embryo_phenotypes = lethality_final[lethality_final["Abnormal Mouse Phenotypes"].str.contains("embryo")]
prefinal_merge = opentargets_merged_no_metals.drop("_merge", axis=1).merge(embryo_phenotypes, on="Gene name", how='outer', indicator=True)
prefinal_filtered = prefinal_merge[prefinal_merge["_merge"] == "left_only"]
print("JAX embryonic lethality filter - unique genes remaining:", prefinal_filtered["Gene name"].unique().shape[0])

#Limit to genes with abnormal inflammation in mice
inflam_df = pd.read_table("data/MGIhdpQuery_markers_20240430_202514_abnormal_inflammatory_response.txt")
inflam_df2 = inflam_df[inflam_df["Organism"]=="mouse"]
inflam_df2["Gene name"] = inflam_df2["Gene Symbol"]
inflam_merged = inflam_df2.merge(mouse_gene_map, on="Gene name")
inflam_final = inflam_merged.merge(mouse_gene_map, on="DB Class Key")
inflam_final["Gene name"] = inflam_final["Gene name_y"]
inflam_merged = inflam_final.merge(prefinal_filtered, on="Gene name")
print("JAX abnormal inflammation filter - unique genes remaining:", inflam_merged["Gene name"].unique().shape[0])

#Bioinformatics filter - only looks at proteins where no other human protein has significant primary sequence homology
#Also filter out proteins that are longer than AlphaFold's max length for 1 protein prediction (and thus split the prediction into multiple files)
bio_df = pd.read_csv("outputs/homology_counts.txt", header=None) #22_with_eval_cut.txt
bio_df["Entry_x"] = bio_df[0].str.split("-").str[1].astype('str')
bio_df["Count"] = bio_df[0].str.split(":").str[-1].astype('int')
bio_df2 = bio_df[bio_df["Count"]==1] #Make sure there are no homologous proteins
bio_df3 =  bio_df2["Entry_x"].value_counts()
bio_df4 =  bio_df3[bio_df3 == 1].reset_index() #Make sure the protein is not longer than AlphaFold's max prediction length
bio_df4.columns = ["Entry_x", "Count2"]
bio_df5 = bio_df4.merge(bio_df2, on="Entry_x")
final_filtered = inflam_merged.merge(bio_df5, on="Entry_x")
print("Bioinformatics filter - unique genes remaining:", final_filtered["Gene name"].unique().shape[0])

# Save final filtered 
final_filtered =  final_filtered.copy()
final_filtered["Entry"] = final_filtered["Entry_x"]
final_filtered[["Gene name", "Entry"]].drop_duplicates().to_csv("outputs/stage1_filtering_gene_candidates.csv", index=None)