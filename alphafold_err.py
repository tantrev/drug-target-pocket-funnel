import Bio.PDB
import pandas as pd
import os

def getGoods(inputPath):
    namey = os.path.splitext(os.path.split(inputPath)[1])[0]
    p = Bio.PDB.PDBParser()
    structure = p.get_structure(namey, inputPath)
    
    listy = list()
    for thing in structure.get_residues():
        res_id = thing.get_id()[1]
        for sub_thing in thing.get_atoms():
            val = sub_thing.get_bfactor()
        listy.append((res_id,val))
    
    df = pd.DataFrame(listy)
    
    goods = df[df[1]>=90]
    return goods