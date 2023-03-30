import pandas as pd
inport numpy as np

from Bio import Seql0
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import os

#rojan's class definition

class OASDBDesc:

    def __init__(self):
        pass

    def read_data(self, rawdata_dir):
        "Gather gz files from the directory and extract these files"

        paired files = [os.path.join(rawdata_dir, f) for f in os.listdir(rawdata_dir) if f.endswith(".gz")]
        t_cols = ['v_call_heavy', 'j_call heavy’, 'v_call_light', 'j_call light', 'sequence_alignment_aa_Tight', 'sequence_alignment aa_heavy','ANARCT_sTatus_Light', 'ANARCI_status_heavy']
        df_seqs = pd.DataFrame()
        for paired_file in paired files:
              df = pd.read_csv(paired_file, compression = ‘gzip’, sep = ', ', skiprows=1)
              df_seqs = pd.concat([df seqs, df[t_cols]], ignore_index=True)
        return df_seqs.copy()

                  
    def perform_random_sample(self, df_seqs, num_iter, n_sample):
        "here we take multiple columns that are random variables"

              "Here two types of data exist - categorical and discrete."
              "we choose one categorical data - v_call heavy"
              "one discrete data - length of sequence alignment_aa_heavy sequence"

              "Since we iterate multiple time to see the sampling is robust,"
              "we create dataframe to store distribution of each sample."

                     

        df_v_heavy = pd.DataFrame()
        df_vh_len = pd.DataFrame()

        "Have a look, if you have time about lambda functions mimic functional programming”""
        df_seqs["VH_Len"] = df_seqs["sequence_alignment_aa_heavy"].apply(lambda row: len(row))

        for i in range(num_iter):
              df_sub_seqs = df_seqs.sample(n_sample)

              "v_call_heavy"
              df_temp = df_sub_seqs[["v_call_heavy"]]
              df temp["iter"]=i
              if df_v_heavy.empty:
                  df_v_heavy = df_temp
              else:
                  df_v_heavy = pd.concat([df_v_heavy, df_temp], ignore_index = True)

                                   

              “heavy chain length"
              df_temp = df_sub_seqs[["VH_Len"]]
              df temp["iter"] = i
              if df_vh_len.empty:
                  df_vh_len= df_temp
              else:
                  df_vh_len = pd.concat([df_vh_len, df_temp], ignore_index = True)

        return df_v_heavy, df_vh_len
   
    def pc_embedding(self, df_sub_sample, seq_col, annotate_col):
        df_selected_sample = df_sub_sample[[annotate col, seq_col]]
        df pc_encode = pd.DataFrame([ProteinAnalysis(i).count_amino_acids() for i in df_selected_sample[seq_col]])
        return df_pc_encode.copy()

    def extended_pc_embedding(self, df_pc_encode)
        df_rsd_mdata = pd. read_csv("rsd_mdata.csv")
        df rsd_mdata.set_index("Aminoacids", inplace=True)

        df_temp = pd.DataFrame()
          # df_temp1 = df_pc_encode[df_rsd_ndata. index]#df_rsd_ndata["PC"].T
          # df temp["PC"] = df_templ.mean(axis=1)

          # df_temp1 = df_pcllencode[df_rsd_ndata. index]*df_rsd_ndata["NC"].T
        # df temp["NC"]= df_templ.mean(axis=1)

                            

          df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["HS"].T
          df temp["HS"] = df_temp1.mean(axis=1)

          df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["pI"].T
          df temp["pI"] = df_temp1.mean(axis=1)

          # df_temp1 = df_pc_encode[df_rsd_ndata. index]*df_rsd_ndata[ "num_atons"].T
          # df temp["num_atoms"] = df_Temp1.mean(axis=1)

          df_temp1 = df_pc_encode[df_rsd_mdata. index]*df_rsd_mdata["hbondDA"].T
          df_temp["hbondDA™] = df_temp1.mean(axis=1)
          return df_temp.join(df_pc_encode, how="inner").copy()

   def pca_analysis(self, df_pc_encode, df_meta, annotate_col)
          #Standard scale
          oscale = StandardScaler()
          # Use fit and transform method
          oscale.fit(df_pc_encode.values)
          encode_scale_data= oscale.transform(df_pc_encode.values)
                  
          #PCA analysis
          opca = PCA(n_components=2)
          opca.fit(encode_scale_data)
          x = opca.transform(encode_scale_data)
          df_pcs = pd.DataFrame(x, columns = ["PC1", "PC2"]

          "Merge PCs with annotation data""”

          df_pcs_meta = df_pcs.join(df_meta, how="inner")
          df pcs_meta[ "newcol"] = df_pcs_meta[annotate_col].apply(lanbda row: row.split("-")[0] \
                                                                                      .split('S')[0]\
                                                                                      .split('D')[0]\
                                                                                      .split('*')[0])
                 

          df_pcs_meta = df_pcs_meta.sort_values("newcol”)
          return df_pcs_meta
                              
                  
                
