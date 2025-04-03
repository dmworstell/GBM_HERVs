#!/usr/bin/env python3
def weighted_random_choice(df):
    import numpy as np
    # Split the dataframe into two: one where exon_number equals 0 and one where it doesn't
    df_zero = df[df['exon_number'] == 0]
    df_non_zero = df[df['exon_number'] != 0]

    exonProb = 1
    
    # If both dataframes are not empty, choose with exonProb probability from df_non_zero and intronProb from df_zero
    if not df_zero.empty and not df_non_zero.empty:
        #return df_non_zero.sample(n=1, weights=np.where(df_non_zero['exon_number']!=0, 0.8, 0.2)).iloc[0] if np.random.rand() < 0.8 else df_zero.sample(n=1).iloc[0]
        return df_non_zero.sample(n=1).iloc[0] if np.random.rand() < exonProb else df_zero.sample(n=1).iloc[0]
    # If df_zero is empty, choose from df_non_zero
    elif df_zero.empty:
        return df_non_zero.sample(n=1).iloc[0]
    # If df_non_zero is empty, choose from df_zero
    elif df_non_zero.empty:
        return df_zero.sample(n=1).iloc[0]