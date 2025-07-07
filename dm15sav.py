import pandas as pd

# Define column names based on your dm15.dat file structure
col_names = ['SN_name', 'SN_dm15', 'SN_dm15source', 'SN_datasource', 'SN_z', ...] 

dm15_df = pd.read_csv(r'./dm15.dat', 
                      names=col_names, 
                      delimiter=',')

# You can now easily access data, e.g.:
# print(dm15_df)
