import pandas as pd
import numpy as np
import glob
import os

from alphafold_err import getGoods

def getPercent(input_x2c,input_df_goods):
    percentList = list()
    for index, row in input_x2c.iterrows():  
        df_site = pd.DataFrame(pd.DataFrame(row["residues"].split(":")[1].strip().split(","))[0].astype('str').str.strip().astype('int')) #s_sitemap_
        truth = df_site.merge(input_df_goods,on=0,how='outer',indicator=True)
        counts = truth["_merge"].value_counts()
        true_denom = counts.sum()-counts[counts.index=="right_only"].iloc[0]
        percent = counts[counts.index=="both"].iloc[0]/true_denom
        minny = df_site[0].min()
        maxxy = df_site[0].max()
        checky = pd.DataFrame(np.arange(minny,maxxy+1))
        checky2 = checky.merge(input_df_goods,indicator=True,how='left')
        all_good_percent = (checky2["_merge"]=="both").sum()/checky2.shape[0]
        left_over = checky2.shape[0]-(checky2["_merge"]=="both").sum()
        percentList.append((percent,all_good_percent,left_over,checky2.shape[0]))
    return percentList

sitemap_file = "outputs/sitemap_results2.csv"
sitemap_base_dir = "pdb_processed/"

df = pd.read_csv(sitemap_file)
cols = ["Entry",]+df.columns[+df.columns.str.contains("sitemap")].tolist() #s_m_entry_name
df["Entry"] =  df["Entry Name"].str.split("-").str[1].astype('str')

subList = list()
for i, file in enumerate(glob.glob(sitemap_base_dir + "/*.pdb")):
    try:
        sub_entry =  os.path.split(file)[1].split("-")[1]
        sub_df = df[df["Entry"]==sub_entry]
        pathy  = file
        goods = getGoods(pathy)
        percents = pd.DataFrame(getPercent(sub_df,goods))
        sub_df["Percent High Quality"] = percents[0].values
        sub_df["Entire Region Percent High Quality"] = percents[1].values
        sub_df["Number Non-High Quality Residues in Entire Region"] = percents[2].values
        sub_df["Total Number of Residues in Entire Region"] = percents[3].values
        subList.append(sub_df)
        print("SUCESSS:",file,i)
    except:
        print("FAILED:",file,i)

end_df = pd.concat(subList)
end_df.to_csv("outputs/sitemap_results_with_quality_metrics.csv", index=None)