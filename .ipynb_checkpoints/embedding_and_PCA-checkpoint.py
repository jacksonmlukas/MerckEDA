import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScalar

import os

class OASDBDesc:
    
    def __init__(self):
        pass
    
    def pc_embedding(self, df_sub_sample, seq_col, annotate_col):
        df_selected_sampled = df_sub_sampled[[annotate_col, seq_col]]
        df_pc_encode = pd.DataFrame([ProteinAnalysis(i).count_amino_acids() for i in df_selected_sample[seq_col]])
        return df_pc_encode.copy()
    
    def extended_pc_embedding(self, df_pc_encode):
        df_rsd_mdata = pd.read_csv("residue_dict_copy.csv")
        df_rsd_mdata.set_index("Aminoacids", inplace=True)
        
        df_temp = pd.DataFrame()
        
        df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["HS"].T
        df_temp["HS"] = df_temp1.mean(axis = 1)
        
        df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["pI"].T
        df_temp["pI"] = df_temp1.mean(axis=1)
        
        df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["hbondDA"].T
        df_temp["hbondDA"] = df_temp1.mean(axis = 1)
        return df_temp.join(df_pc_encode, how="inner").copy()
    
    def pca_analysis(self, df_pc_encode, df_meta, annotate_col):
        
        oscale = StandardScaler()
        oscale.fit(df_pc_encode.values)
        encode_scale_data = oscale.transform(df_pc_encode.values)
        
        opca = PCA(n_components=2)
        opca.fit(encode_scale_data)
        x = opca.transform(encode_scale_data)
        df_pcs = pd.DataFrame(x, columns = ["PC1", "PC2"])
        
        df_pcs_meta = df_pcs.join(df_mera, how="inner")
        df_pcs_meta["newcol"] = df_pcs_meta[annotate_col].apply(lambda row: row.split("-")[0] \
                                                                                .split('S')[0]\
                                                                                .split('D')[0]\
                                                                                .split("*")[0])
        df_pcs_meta = df_pcs_meta.sort_values("newcol")
        return df_pcs_meta