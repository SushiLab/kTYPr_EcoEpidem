import numpy as np
import pandas as pd

def group_genomes_by_ani(matrix_file, ani_threshold):
    with open(matrix_file, 'r') as f:
        lines = f.readlines()[1:]  # Ignore the first line
    
    genomes = []
    ani_matrix = []

    for i, line in enumerate(lines):
        parts = line.strip().split()
        genomes.append(parts[0].split('/')[-2])  # First column is genome path
        ani_values = [float(x) for x in parts[1:]]  # Convert ANI values to float
        ani_matrix.append(ani_values)

    # Convert to a square matrix by filling missing values with NaN
    max_len = len(ani_matrix)
    ani_matrix = [row + [np.nan] * (max_len - len(row)) for row in ani_matrix]
    ani_matrix = np.array(ani_matrix)

    representatives = []
    genome_groups = {}

    for i, genome in enumerate(genomes):
        if not representatives:
            representatives.append(genome)
            genome_groups[genome] = [genome]
        else:
            best_rep = None
            best_avg_ani = 0

            for rep in representatives:
                rep_index = genomes.index(rep)
                if i < rep_index:
                    ani_value = ani_matrix[rep_index, i]
                else:
                    ani_value = ani_matrix[i, rep_index]

                if not np.isnan(ani_value) and ani_value >= ani_threshold:
                    avg_ani = np.nanmean([ani_matrix[max(i, j), min(i, j)] for j in range(len(genomes)) if genomes[j] in genome_groups[rep]])
                    if avg_ani > best_avg_ani:
                        best_avg_ani = avg_ani
                        best_rep = rep

            if best_rep:
                genome_groups[best_rep].append(genome)
            else:
                representatives.append(genome)
                genome_groups[genome] = [genome]

    return genome_groups

def reverse_and_merge_dicts(dict98, dict99):
    # Reverse each dictionary
    reversed_98 = {genome: rep for rep, genomes in dict98.items() for genome in genomes}
    reversed_99 = {genome: rep for rep, genomes in dict99.items() for genome in genomes}
    
    # Get all unique genomes as a sorted list
    all_genomes = sorted(list(reversed_98.keys()))
    
    # Create a DataFrame
    df = pd.DataFrame({'genome': all_genomes})
    df['ANI98'] = df['genome'].map(reversed_98)
    df['ANI99'] = df['genome'].map(reversed_99)
    
    return df

file_path = 'skani_matrix.txt'   # Output from ./skani_derep.sh
anis98 = group_genomes_by_ani(file_path, 98)
anis99 = group_genomes_by_ani(file_path, 99)
df = reverse_and_merge_dicts(anis98, anis99)
df.to_csv('ani98_ani99_representatives.tsv', index=False, sep='\t')