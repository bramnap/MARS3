import pandas as pd
from . import agora_checking

def taxonomic_distribution(total_reads, absent, present, agora2_level_set, df, level, levels):
    
    df_present_species = df.loc[df['Species'] != pd.NA] #present rows in general dataframe
    levels_copy = levels.copy()
    levels_copy.remove(level)
    levels_copy.append("Kingdom") #for dropping
    df_present_species = df_present_species.drop(columns=levels_copy).groupby(level).sum()
    df_absent_agora_species, df_present_agora_species = agora_checking.agora_checking(df_present_species, agora2_level_set)
    
    return df_present_species, df_present_agora_species, df_absent_agora_species, absent/total_reads