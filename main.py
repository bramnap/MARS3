from lib import preprocessing, agora_checking, taxonomic_distribution, pipeline, species_genus_association, general_stats

import os
import json

def main(*args, relative=False, **kwargs):
    with open('mars.json', 'r') as fp:
        agora2 = json.load(fp)
    agora2_phyla = set(agora2["Phylum"])
    agora2_classes = set(agora2["Class"])
    agora2_orders = set(agora2["Order"])
    agora2_families = set(agora2["Family"])
    agora2_genera = set(agora2["Genus"])
    agora2_species = set(agora2["Species"])
    agora2_strains = set(agora2["Strain"])
    
    agora2_level_sets = [agora2_phyla, agora2_classes, agora2_orders, agora2_families, agora2_genera, agora2_species, agora2_strains]
    
    df_levels = list(preprocessing.preprocessing(**kwargs, relative=relative))
    # Total reads
    df, kingdom_df, phylum_df, class_df, order_df, family_df, genus_df, species_df, strain_df = df_levels[0], df_levels[1], df_levels[2], df_levels[3], df_levels[4], df_levels[5], df_levels[6], df_levels[7], df_levels[8]

    # phylum_df total phylum reads
    
    # associated_phylum_species, associated_phylum_genus = species_genus_association.association(df, levels, "Phylum")
    # absent_phylum, present_phylum = agora2_checking(phylum_df, agora2_phyla) 
    # associated_phylum_with_species in agora

    present_df_levels = []
    absent_df_levels = []
    for df_level, agora2_level_set in zip(df_levels[2:], agora2_level_sets):
        df_absent, df_present = agora_checking.agora_checking(df_level, agora2_level_set)
        present_df_levels.append(df_present)
        absent_df_levels.append(df_absent)
        
#     present_phylum_df, present_class_df, present_order_df, present_family_df, present_genus_df, present_species_df, present_strain_df = present_df_levels[0], present_df_levels[1], present_df_levels[2], present_df_levels[3], present_df_levels[4], present_df_levels[5], present_df_levels[6]
#     absent_phylum_df, absent_class_df, absent_order_df, absent_family_df, absent_genus_df, absent_species_df, absent_strain_df = absent_df_levels[0], absent_df_levels[1], absent_df_levels[2], absent_df_levels[3], absent_df_levels[4], absent_df_levels[5], absent_df_levels[6]
    present_genus_df, present_species_df = present_df_levels[4], present_df_levels[5]
    
    # #construct coverage files here
    # levels_omitting_kingdom = levels[1:].copy()
    # total_reads = df.drop(columns=levels).sum()
    
    # present_dataframes = {}
    # absent_dataframes = {}
    # for i, (absent_df, present_df, level, agora2_level_set) in enumerate(zip(absent_df_levels, present_df_levels, levels_omitting_kingdom, agora2_level_sets)):
    #     df_present_species, df_present_agora_species, df_absent_species, absent_relative = taxonomic_distribution.taxonomic_distribution(total_reads, absent_df, present_df, agora2_level_set, df, level, levels_omitting_kingdom)
    #     present_dataframes[level.lower()] = [present_df_levels[i], df_present_species, df_present_agora_species]
    #     absent_dataframes[level.lower()] = [absent_df_levels[i], df_absent_species, absent_relative]
    
    levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

    if not os.path.isdir("MARS_output"):
        os.mkdir("MARS_output")

    #save requested files
    for arg in args:
        try:
            arg = arg.lower()
            if arg == "class":
                #save
                #print(present_dataframes[arg])
                pass
            elif arg == "order":
                #save
                pass
            elif arg == "family":
                #save
                pass
            elif arg == "strain":
                #save
                pass
            else:
                print(f"\"{arg}\" did not match any of the optional taxonomic levels. Please check spelling")
        except SyntaxError:
            print("A syntax error was found in your arguments. Please check that you inputted a string.")
            pass


    # total_df, associated_species, associated_genus, absent, present, associated_species_agora2, associated_genus_agora2 = pipeline.pipeline(df, total_df, levels, level, agora2_level_set, agora2_species, agora2_genera)
    species_phylum_list, genus_phylum_list = pipeline.pipeline(df, phylum_df, levels, "Phylum", agora2_phyla, agora2_species, agora2_genera)

    # for taxa in ['phylum', 'genus', 'species']:
    #     for i, name in enumerate(["agora_checked", "total_with_species", "agora2"]):
    #         present_dataframes[taxa][i].to_csv(f'MARS_output/{name}_{taxa}_present.csv')
    #     for i, name in enumerate(["agora_checked", "agora2", "relative"]):
    #         absent_dataframes[taxa][i].to_csv(f'MARS_output/{name}_{taxa}_absent.csv')

    # Normalise
    # TODO: might want to put this in a separate module?

    total_species_reads = present_species_df.sum()
    total_genus_reads = present_genus_df.sum()

    agora_species_normed = present_species_df.loc[:].div(total_species_reads)
    agora_genus_normed = present_genus_df.loc[:].div(total_genus_reads)

    # Save these dfs
    agora_species_normed[agora_species_normed < 1e-5] = 0
    agora_genus_normed[agora_genus_normed < 1e-5] = 0

    # Renormalize

    total_species_rel_abund = agora_species_normed.sum()
    total_genus_rel_abund = agora_genus_normed.sum()


    # Get stats on the species

    x = general_stats.general_stats(df, species_phylum_list, present_species_df)

    # Get stats on the genus

    return present_genus_df, present_species_df

if __name__ == "__main__":

    genus, species = main(taxonomy_table=r"C:\Users\MSPG\Desktop\Mars_test\taxonomyWoL.tsv",
                          feature_table=r"C:\Users\MSPG\Desktop\Mars_test\feature-tableWoLgenome.txt")

