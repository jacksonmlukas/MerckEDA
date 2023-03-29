import pandas as pd

import ablang
import numpy as np
import os
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
import matplotlib
import matplotlib.pyplot as plt

#rojan's class definition

class OASDBDesc:
    
    def __init__(self):
        pass

    def read_data(self, rawdata_dir):
        "Gather gz files from the directory and extract these files"
    
        paired_files = [os.path.join(rawdata_dir, f) for f in os.listdir(rawdata_dir) if f.endswith(".gz")] 
        t_cols = ['v_call_heavy', 'd_call_heavy', 'j_call_heavy', 'sequence_alignment_aa_light', 
                  'sequence_alignment_aa_heavy', 'ANARCI_status_light', 'ANARCI_status_heavy']

        df_seqs = pd.DataFrame()
        for paired_file in paired_files:
            df = pd.read_csv(paired_file, compression = 'gzip', sep=',', skiprows=1)
            df_seqs = pd.concat([df_seqs, df[t_cols]], ignore_index=True)
        return df_seqs.copy()

    def encode_seq(self, df, column):
        #function to encode sequences
        
        #5 main types of protein encoding methods: binary encoding, 
        #physiochemical properties encoding, evolution-based encoding, structure-based encoding, 
        #and machine-learning encoding.
        
        #ablang
        
        #heavy sequence encoding
        heavy_ablang = ablang.pretrained("heavy")
        heavy_ablang.freeze()
        
        seqs_heavy = df.loc[1:30, 'sequence_alignment_aa_heavy']

        seqcodings_heavy = heavy_ablang(seqs_heavy, mode='seqcoding')
        print("-"*100)
        print("The output shape of the heavy seq-codings:", seqcodings_heavy.shape)
        print("-"*100)

        print(seqcodings_heavy)
        
        #light sequence encoding
        light_ablang = ablang.pretrained("light")
        light_ablang.freeze()
        
        seqs_light = df.loc[1:30, 'sequence_alignment_aa_light']

        seqcodings_light = light_ablang(seqs_light, mode='seqcoding')
        print("-"*100)
        print("The output shape of the light seq-codings:", seqcodings_light.shape)
        print("-"*100)

        print(seqcodings_light)
        
    
    #cheryl's one-hot encoding 
        
    def one_hot_encode_seq(self, df, column):
    #Output a df with a specific columns that want to get dummies in
    
        #label_encode
        le = LabelEncoder()
        le.fit(df[column])
        integer_encoded_letters_arry = le.transform(df[column])

        #append
        integer_encoded_letters_series = pd.Series(integer_encoded_letters_arry)
        df['integer_encoded_letters'] = integer_encoded_letters_series

        #one hot encode
        df_dummies = pd.get_dummies(df, prefix = ['integer_encoded_letters'], columns = ['integer_encoded_letters'], drop_first = True)
        #return df_dummies
    
        
    #cheryl's code - physiochemical properties encoding
    
    def physchemvh_gen(self, df, column):
        alph = np.array(sorted('ACDEFGHIKLMNPQRSTVWY'))
        residue_info = pd.read_csv("residue_dict_copy.csv", header = 0, index_col = 0)
        
        res_counts = pd.DataFrame(index = alph)
        df = df.set_index(column)
        for i in df.index:
            characters = pd.Series(list(i))
            res_counts = pd.concat([res_counts, characters.value_counts()], axis = 1, ignore_index = False)
        res_counts.fillna(0, inplace = True)
        res_counts = res_counts.T
        hydrophobicity = []    
        for column in res_counts:
            hydros = []
            for index, row in res_counts.iterrows():
                hydros.append(row[column]*residue_info.loc[column, 'Hydropathy Score'])
            hydrophobicity.append(hydros)
        hydrophobicity = pd.DataFrame(hydrophobicity).T
        hydrophobicity['ave'] = hydrophobicity.sum(axis = 1)/115
        res_counts['Hydro'] = res_counts['A'] +  res_counts['I'] +  res_counts['L']+  res_counts['F']+  res_counts['V']
        res_counts['Amph'] = res_counts['W'] +  res_counts['Y']+  res_counts['M']
        res_counts['Polar'] = res_counts['Q'] +  res_counts['N'] + res_counts['S'] +  res_counts['T'] +  res_counts['C']+  res_counts['M']
        res_counts['Charged'] =  res_counts['R'] +  res_counts['K'] + res_counts['D'] +  res_counts['E'] +  res_counts['H']
        res_counts.reset_index(drop = True, inplace = True)
        physchemvh = pd.concat([res_counts, hydrophobicity['ave']], axis = 1, ignore_index = False)
        
        return "One Hot Encoding: ", df_dummies, "Physiochemical Properties Encoding: ", physchemvh


# In[8]:


# how to use class
#rawdata_dir = "/datasets/merck-files/Merck-files/"
obj = OASDBDesc()
#df_seq = obj.read_data(rawdata_dir)


# In[ ]:
df_seq =  pd.read_csv("/datasets/merck-files/Merck-files/sampled_df.csv")

#df_seq.head() 

#annotate all of these columns


# In[ ]:

#obj.encode_seq(df_seq, "sequence_alignment_aa_heavy")




