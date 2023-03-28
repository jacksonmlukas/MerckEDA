import pandas as pd
inport numpy as np

import os

#rojan's class definition

class OASDBDesc:

    def _init_(self):
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
