# =======================================
########### Input functions #############
# =======================================

def get_samples_for_cond(con):
    
    # return list of all samples for the current condition
    samples = df_cond.loc[df_cond['condition']==con]['sample'].tolist()
    
    return samples

def get_majiq_deltapsi_group_id(comp):

    # Get cond1 and cond2 as a list
    return comp.split("-")

def get_majiq_files_deltapsi(comp):

    # Get all .majiq files for the 2 conditions to compare
    cond_list = get_majiq_deltapsi_group_id(comp)
    cond1, cond2 = cond_list[0], cond_list[1]

    cond1_samples = get_samples_for_cond(cond1)
    cond2_samples = get_samples_for_cond(cond2)

    cond1_majiq = [
        'results/splicing/majiq/build/{id}_Aligned.sortedByCoord.out.primary.majiq'
        for id in cond1_samples
    ]
    cond2_majiq = [
        'results/splicing/majiq/build/{id}_Aligned.sortedByCoord.out.primary.majiq'
        for id in cond2_samples
    ]
    return cond1_majiq, cond2_majiq


def majiq_cond1(wildcards):

    # Return .majiq files for 1st condition only
    cond1_majiq, cond2_majiq = get_majiq_files_deltapsi(wildcards.comp)
    return cond1_majiq

def majiq_cond2(wildcards):

    # Return .majiq files for 2nd condition only
    cond1_majiq, cond2_majiq = get_majiq_files_deltapsi(wildcards.comp)
    return cond2_majiq
