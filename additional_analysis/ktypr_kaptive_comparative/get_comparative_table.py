# Load the different ktypr files

import pandas as pd

ktypr = {'kref_wg': pd.read_csv("./to_report/wg_flank_v2/kref_wg_20241203_ktyps.tsv", sep='\t'),
         'kr10_wg': pd.read_csv("./to_report/wg_flank_v2/kref10_wg_20241203_ktyps.tsv", sep='\t'),
         'binf_wg': pd.read_csv("./tool_comparative/comp_data/ktypr_ktypr/results_ktypr.tsv", sep='\t'), 
         'kref_fl': pd.read_csv("./to_report/wg_flank_v2/kref_flank_20241209_ktyps.tsv", sep='\t'),
         'kr10_fl': pd.read_csv("./to_report/wg_flank_v2/kref10_flank_20241209_ktyps.tsv", sep='\t'),
         'binf_fl': pd.read_csv("./tool_comparative/comp_data/ktypr_ktypr/results_ktypr.tsv", sep='\t'),
         'kapt_wg': pd.read_csv("./tool_comparative/comp_data/kaptive_ktypr/ktypr_wgresults_ktypr.tsv", sep='\t'),                  
         'kapt_fl': pd.read_csv("./tool_comparative/comp_data/kaptive_ktypr/ktypr_wgresults_ktypr.tsv", sep='\t'),                  
         }
ktypr['kapt_fl'].index = [i.replace('ktypr_wg', '') for i in ktypr['kapt_fl'].index]
ktypr['kapt_wg'].index = [i.replace('ktypr_wg', '') for i in ktypr['kapt_wg'].index]

kaptive = {'kref_kp': pd.read_csv("./tool_comparative/comp_data/kaptive_kref.tsv", sep='\t'),
           'kr10_kp': pd.read_csv("./tool_comparative/comp_data/kaptive_kref10.tsv", sep='\t'),
           'binf_kp': pd.read_csv("./tool_comparative/comp_data/kaptive_bioinforef.tsv", sep='\t'),
           'kapt_kp': pd.read_csv("./tool_comparative/comp_data/kaptive_kaptive.tsv", sep='\t')}

kaptive_kaps = {'kref_kp_kdb': pd.read_csv("./tool_comparative/comp_data/kaptive_ktyprDB_kref.tsv", sep='\t'),
                'kr10_kp_kdb': pd.read_csv("./tool_comparative/comp_data/kaptive_ktyprDB_kref10.tsv", sep='\t'),
                'binf_kp_kdb': pd.read_csv("./tool_comparative/comp_data/kaptive_ktyprDB_bioinforef.tsv", sep='\t'),
                'kapt_kp_kdb': pd.read_csv("./tool_comparative/comp_data/kaptive_ktyprDB_kaptive.tsv", sep='\t')}

# Get the best for kaptive
for k, v in kaptive.items():
    for _, row in v.iterrows():
        genome = row['Assembly']
        kaptive_type = row['Best match locus']
        rs[genome]+=[kaptive_type, genome.split('_')[0]]

# Build table
rs = {}

fl26 = pd.read_csv("./comp_data/ktypr26/flresults_ktypr.tsv", sep='\t')
wg26 = pd.read_csv("./comp_data/ktypr26/wgresults_ktypr.tsv", sep='\t')

cols = ['dataset', 'ktypr_wg', 'ktypr_flank', 'kaptive','kaptive_with_KapsDB', 'expected', 'ktypr_wg_2026', 'kytpr_fl_2026']
for ty in ['kref', 'kr10', 'binf', 'kapt']:
    for flank in ['wg', 'fl']:
        run = f"{ty}_{flank}"
        for _, row in ktypr[run].iterrows():
            try:
                genome = row['genome_id']
            except:
                genome = _
            ktypr_type = row['predicted']
            if genome not in rs:
                rs[genome] = [ty, ktypr_type]
            else:
                rs[genome].append(ktypr_type)

for ty in ['kref', 'kr10', 'binf', 'kapt']:
    run = f"{ty}_kp"
    seen_genomes = set()
    for _, row in kaptive[run].iterrows():
        try:
            genome = row['Assembly']
        except:
            genome = _
        kaptive_type = row['Best match locus']
        if genome not in seen_genomes:
            rs[genome]+= [kaptive_type]
            seen_genomes.add(genome)

for ty in ['kref', 'kr10', 'binf', 'kapt']:
    run = f"{ty}_kp_kdb"
    seen_genomes = set()
    for _, row in kaptive_kaps[run].iterrows():
        try:
            genome = row['Assembly']
        except:
            genome = _
        kaptive_type = row['Best match locus']
        if genome not in seen_genomes:
            rs[genome]+= [kaptive_type, genome.split('_')[0]]
            seen_genomes.add(genome)

for flank, df in zip(['wg', 'fl'], [wg26, fl26]):
    for _, row in df.iterrows():
        genome = _.replace('wg', '').replace('fl', '')
        ktypr_type = row['predicted']
        if genome in rs:
            rs[genome].append(ktypr_type)
        else:
            print(genome)

results_df = pd.DataFrame.from_dict(rs, orient='index', columns=cols)
results_df.reset_index(inplace=True)
results_df.rename(columns={'index': 'genome_id'}, inplace=True)

results_df.to_excel("/nfs/home/smiravet/KTYPS_DEV/scratch/tool_comparative/comp_data/comparative_ktypr_kaptive_20250116.xlsx", index=False)