import pandas as pd
import numpy as np

from Bio import SeqIO
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

        paired_files = [os.path.join(rawdata_dir, f) for f in os.listdir(rawdata_dir) if f.endswith(".gz")]
        t_cols = ['v_call_heavy', 'j_call_heavy', 'v_call_light', 'j_call_light', 'sequence_alignment_aa_light', 'sequence_alignment_aa_heavy','ANARCI_status_light', 'ANARCI_status_heavy']
        df_seqs = pd.DataFrame()
        for paired_file in paired_files:
            print(paired_file)
            df = pd.read_csv(paired_file, compression = "gzip", sep = ",", skiprows=1)
            df_seqs = pd.concat([df_seqs, df[t_cols]], ignore_index=True)
        return df_seqs.copy()

                  
    def perform_random_sample(self, df_seqs, num_iter, n_sample):
        #here we take multiple columns that are random variables

              #Here two types of data exist - categorical and discrete."
              #"we choose one categorical data - v_call heavy"
              #"one discrete data - length of sequence alignment_aa_heavy sequence"

             # "Since we iterate multiple time to see the sampling is robust,"
            #  "we create dataframe to store distribution of each sample."

                     

        df_v_heavy = pd.DataFrame()
        df_vh_len = pd.DataFrame()

        #Have a look, if you have time about lambda functions mimic functional programming”""
        df_seqs["VH_Len"] = df_seqs["sequence_alignment_aa_heavy"].apply(lambda row: len(row))

        for i in range(num_iter):
            df_sub_seqs = df_seqs.sample(n_sample)

              #"v_call_heavy"
            df_temp = df_sub_seqs[["v_call_heavy"]]
            df_temp["iter"]=i
            if df_v_heavy.empty:
                df_v_heavy = df_temp
            else:
                df_v_heavy = pd.concat([df_v_heavy, df_temp], ignore_index = True)

                                   

              #“heavy chain length"
        df_temp = df_sub_seqs[["VH_Len"]]
        df_temp["iter"] = i
        if df_vh_len.empty:
            df_vh_len= df_temp
        else:
            df_vh_len = pd.concat([df_vh_len, df_temp], ignore_index = True)

        return df_v_heavy, df_vh_len
   
    def pc_embedding(self, df_sub_sample, seq_col, annotate_col):
        df_selected_sample = df_sub_sample[[annotate_col, seq_col]]
        df_pc_encode = pd.DataFrame([ProteinAnalysis(i).count_amino_acids() for i in df_selected_sample[seq_col]])
        return df_pc_encode.copy()

    def extended_pc_embedding(self, df_pc_encode):
        df_rsd_mdata = pd.read_csv("rsd_mdata.csv")
        df_rsd_mdata.set_index("Aminoacids", inplace=True)

        df_temp = pd.DataFrame()
          # df_temp1 = df_pc_encode[df_rsd_ndata. index]#df_rsd_ndata["PC"].T
          # df temp["PC"] = df_templ.mean(axis=1)

          # df_temp1 = df_pcllencode[df_rsd_ndata. index]*df_rsd_ndata["NC"].T
        # df temp["NC"]= df_templ.mean(axis=1)

                            

        df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["HS"].T
        df_temp["HS"] = df_temp1.mean(axis=1)

        df_temp1 = df_pc_encode[df_rsd_mdata.index]*df_rsd_mdata["pI"].T
        df_temp["pI"] = df_temp1.mean(axis=1)

          # df_temp1 = df_pc_encode[df_rsd_ndata. index]*df_rsd_ndata[ "num_atons"].T
          # df temp["num_atoms"] = df_Temp1.mean(axis=1)

        df_temp1 = df_pc_encode[df_rsd_mdata. index]*df_rsd_mdata["hbondDA"].T
        df_temp["hbondDA"] = df_temp1.mean(axis=1)
        return df_temp.join(df_pc_encode, how="inner").copy()

    def pca_analysis(self, df_pc_encode, df_meta, annotate_col):
          #Standard scale
        oscale = StandardScaler()
          # Use fit and transform method
        oscale.fit(df_pc_encode.values)
        encode_scale_data= oscale.transform(df_pc_encode.values)
                  
          #PCA analysis
        opca = PCA(n_components=2)
        opca.fit(encode_scale_data)          
        x = opca.transform(encode_scale_data)
        df_pcs = pd.DataFrame(x, columns = ["PC1", "PC2"])

          #"Merge PCs with annotation data""”

        df_pcs_meta = df_pcs.join(df_meta, how="inner")
        df_pcs_meta["newcol"] = df_pcs_meta[annotate_col].apply(lambda row: row.split("-")[0] \
                                                                                      .split('S')[0]\
                                                                                      .split('D')[0]\
                                                                                      .split('*')[0])
                 

        df_pcs_meta = df_pcs_meta.sort_values("newcol")
        return df_pcs_meta
    
    #one hot encode
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
    return df_dummies
    
    #physchemvh_gen
    alph = np.array(sorted('ACDEFGHIKLMNPQRSTVWY'))
    residue_info = pd.read_csv("residue_dict_copy.csv", header = 0, index_col = 0)
    def physchemvh_gen(self, df, column):
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
        hydrophobicity['ave'] = hydrophobicity.mean()
        res_counts['Hydro'] = res_counts['A'] +  res_counts['I'] +  res_counts['L']+  res_counts['F']+  res_counts['V']
        res_counts['Amph'] = res_counts['W'] +  res_counts['Y']+  res_counts['M']
        res_counts['Polar'] = res_counts['Q'] +  res_counts['N'] + res_counts['S'] +  res_counts['T'] +  res_counts['C']+  res_counts['M']
        res_counts['Charged'] =  res_counts['R'] +  res_counts['K'] + res_counts['D'] +  res_counts['E'] +  res_counts['H']
        res_counts.reset_index(drop = True, inplace = True)
        physchemvh = pd.concat([res_counts, hydrophobicity['ave']], axis = 1, ignore_index = False)
        return physchemvh
    
    #find the best cluster
    def best_num_cluster(X, elbow = True, silhouette = True):
    #Elbow Method
    if elbow == True:
        Sum_of_squared_distances = []
        K = range(1,len(X) + 1)
        for num_clusters in K :
            kmeans = KMeans(n_clusters=num_clusters, random_state = 48)
            kmeans.fit(X)
            Sum_of_squared_distances.append(kmeans.inertia_)
        plt.plot(K,Sum_of_squared_distances,'bx-')
        plt.xlabel('Values of K') 
        plt.ylabel('Sum of squared distances/Inertia') 
        plt.title('Elbow Method For Optimal k')
        plt.show()
        
    if silhouette == True:
        #Silhouette
        range_n_clusters = range(2, len(X))
        silhouette_avg = []
        for num_clusters in range_n_clusters:
             # initialise kmeans
            kmeans = KMeans(n_clusters=num_clusters, random_state = 48)
            kmeans.fit(X)
            cluster_labels = kmeans.labels_

            # silhouette score
            silhouette_avg.append(silhouette_score(X, cluster_labels))
        plt.plot(range_n_clusters,silhouette_avg,'bx-')
        plt.xlabel('Values of K') 
        plt.ylabel('Silhouette score') 
        plt.title('Silhouette analysis For Optimal k')
        plt.show()
                              
                  
                
