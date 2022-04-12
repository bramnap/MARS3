
def normalise_and_cut(present_df):

    total_reads = present_df.sum()
    agora_normed = present_df.loc[:].div(total_reads)

    agora_normed_cut = agora_normed.copy()
    agora_normed_cut[agora_normed_cut < 1e-5] = 0

    # Renormalize
    total_rel_abund = agora_normed_cut.sum()
    agora_renormed = agora_normed_cut.loc[:].div(total_rel_abund)

    return agora_normed, agora_normed_cut, agora_renormed
