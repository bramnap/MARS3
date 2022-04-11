import pandas as pd


def general_stats(initial_df, list_phylum_df, agora2_df):

    total_phylum, associated_phylum, associated_agora_phylum = list_phylum_df

    sum_initial_df = sum_rename_sort(initial_df.iloc[:,8:], 'Total Reads')

    sum_associated_reads = sum_rename_sort(associated_phylum,'Associated Reads')
    final_df = pd.merge(sum_initial_df, sum_associated_reads, left_index=True, right_index=True)

    sum_associated_agora_reads = sum_rename_sort(associated_agora_phylum,'AGORA2 associated reads')
    final_df = pd.merge(final_df, sum_associated_agora_reads, left_index=True, right_index=True)

    reads_cov = sum_associated_reads/sum_initial_df
    reads_cov = reads_cov.rename('% Reads that have species')
    final_df = pd.merge(final_df, reads_cov, left_index=True, right_index=True)

    reads_agora_cov = sum_associated_agora_reads / sum_initial_df
    reads_agora_cov = reads_agora_cov.rename('% Agora2 total reads coverage')
    final_df = pd.merge(final_df, reads_agora_cov, left_index=True, right_index=True)

    reads_agora_cov = sum_associated_agora_reads / sum_associated_reads
    reads_agora_cov = reads_agora_cov.rename('% Agora2 associated reads coverage')
    final_df = pd.merge(final_df, reads_agora_cov, left_index=True, right_index=True)



    print(final_df)
    return


def sum_rename_sort(df, name):
    df_sum = df.sum(axis=0)
    df_sumname = df_sum.rename(name)
    df_f = df_sumname.sort_index()

    return df_f
