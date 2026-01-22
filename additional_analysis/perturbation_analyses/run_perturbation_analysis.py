import os
import glob
import random
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import collections
import pyhmmer


GBK_DIR = "./kTYPr/ktypr/data/reference_clusters"   # Path to GenBank files from reference clusters zip in kTYPr repo

rows = []

for fname in sorted(os.listdir(GBK_DIR)):
    if not fname.endswith((".gbk", ".gbff")):
        continue

    path = os.path.join(GBK_DIR, fname)

    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            if feature.location is None:
                continue

            # Extract nucleotide sequence using coordinates
            nt_seq = feature.extract(record.seq).upper()
            nt_seq = str(nt_seq)

            if len(nt_seq) == 0:
                continue

            gc = 100 * (nt_seq.count("G") + nt_seq.count("C")) / len(nt_seq)

            gene = feature.qualifiers.get("gene", ["NA"])[0]
            hmm  = feature.qualifiers.get("hmm_id", ["NA"])[0]

            rows.append({
                "genbank": fname,
                "gene": gene,
                "hmm_id": hmm,
                "gc": gc
            })

df = pd.DataFrame(rows)

#########################
#### GC exploration plots

import matplotlib.pyplot as plt
import numpy as np
import math

# --- Parameters
genbanks = df["genbank"].unique()
per_row = 20
n_rows = math.ceil(len(genbanks) / per_row)

fig, axes = plt.subplots(n_rows, 1, figsize=(max(8, 0.6*per_row), 4*n_rows), sharey=True)

if n_rows == 1:
    axes = [axes]  # make it iterable

for row_idx in range(n_rows):
    ax = axes[row_idx]
    start = row_idx * per_row
    end = min((row_idx+1) * per_row, len(genbanks))
    gb_subset = genbanks[start:end]

    for i, gb in enumerate(gb_subset):
        sub = df[df["genbank"] == gb]
        x = np.random.normal(i, 0.08, size=len(sub))  # jitter
        ax.scatter(x, sub["gc"], alpha=0.7)

    ax.set_xticks(range(len(gb_subset)))
    ax.set_xticklabels(gb_subset, rotation=90)
    ax.set_ylabel("GC content (%)")
    ax.set_xlabel("GenBank file")

plt.suptitle("Per-gene GC content across reference clusters", y=1.02)
plt.tight_layout()
plt.show()
fig.savefig("gc_content_reference_clusters.png", dpi=300, bbox_inches='tight')


#########################
#### Perturbation (mutation and trimming) analysis code

# --- Step 0: protein alphabet (needed for retrieve_hits)
protein_alphabet = pyhmmer.easel.Alphabet.amino()

# --- Step 1: extract proteins from genbank files
def extract_proteins_from_gbks(gbk_files):
    proteins = {}
    for gbk in gbk_files:
        for record in SeqIO.parse(gbk, "genbank"):
            for feat in record.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    locus = feat.qualifiers.get("hmm_id", [None])[0]
                    seq = feat.qualifiers["translation"][0]
                    if locus and locus not in proteins:
                        proteins[locus] = seq
    return proteins

# --- Step 2: generate degraded/mutated sequences
AA = list("ACDEFGHIKLMNPQRSTVWY")

def degrade_n(seq, frac): n = int(len(seq)*frac); return seq[n:]
def degrade_c(seq, frac): n = int(len(seq)*frac); return seq[:-n] if n>0 else seq
def degrade_both(seq, frac): n = int(len(seq)*frac); return seq[n:-n] if n>0 else seq
def mutate(seq, identity):
    seq = list(seq)
    n_mut = int(len(seq)*(1-identity))
    for i in random.sample(range(len(seq)), n_mut):
        seq[i] = random.choice([a for a in AA if a != seq[i]])
    return "".join(seq)

def write_fasta(seqs_dict, path):
    records = [SeqRecord(Seq(seq), id=locus, description="") for locus, seq in seqs_dict.items()]
    SeqIO.write(records, path, "fasta")

# --- Step 3: generate variants and run HMMs
def run_variants(proteins, hmms_path, tmp_dir="./tmp_seqs"):
    os.makedirs(tmp_dir, exist_ok=True)
    results_all = []

    fracs = [0.1,0.2,0.3,0.4]
    identities = [i/10 for i in range(10,-1,-1)]

    for locus, seq in proteins.items():
        variants = {}

        # N-term, C-term, Both
        for f in fracs:
            variants[f"N_{int(f*100)}"] = degrade_n(seq, f)
            variants[f"C_{int(f*100)}"] = degrade_c(seq, f)
            variants[f"Both_{int(f*100)}"] = degrade_both(seq, f)

        # Mutations
        for ident in identities:
            variants[f"Mut_{int(ident*100)}"] = mutate(seq, ident)

        # Write variants to fasta
        fasta_path = os.path.join(tmp_dir, f"{locus}.fasta")
        write_fasta(variants, fasta_path)

        # Run HMMs
        df_hits = retrieve_hits(fasta_path, hmms_path)
        df_hits["locus"] = locus
        results_all.append(df_hits)

    # Combine
    return pd.concat(results_all, ignore_index=True)

# --- Step 4: Plot bitscore impact
def plot_bitscore(df, locus):
    sub = df[df.locus==locus]
    fig, ax = plt.subplots(figsize=(8,5))

    for mode in sub.query.unique():
        d = sub[sub.query==mode]
        ax.plot(d.subject, d.bitscore, marker="o", label=mode)

    ax.set_xlabel("Variant")
    ax.set_ylabel("Bitscore")
    ax.set_title(locus)
    ax.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

# Retrieve hits

def retrieve_hits(seqs_path, hmms_path, 
                  fields=["query", "subject", "bitscore", "evalue"],  max_aa_size=2000):
    """
    Given a sequence path <seqs_path> in fasta format and a directory with hmms <hmms_path>,
    returns a dataframe with the hits containing <fields>
    """
    # Load cluster proteins
    with pyhmmer.easel.SequenceFile(seqs_path, digital=True, alphabet=protein_alphabet) as seqs_file:
        proteins = [seq for seq in list(seqs_file) if len(seq) <= max_aa_size] 

    # Load HMMs
    hmms = []
    for fil in glob.glob(hmms_path):
        with pyhmmer.plan7.HMMFile(fil) as hmm_file:
            hmms.append(hmm_file.read())

    # Run HMMs
    Result = collections.namedtuple("Result", fields)

    results = []
    for hits in pyhmmer.hmmsearch(hmms, proteins, E=1):
        cog = hits.query_name.decode()
        for hit in hits:
            if hit.included:
                results.append(Result(hit.name.decode(), cog, hit.score, hit.evalue))

    # Results --> df
    hits_df = {}
    c = 0
    for i in results:
        hits_df[c] = list(i)
        c += 1
    hits_df = pd.DataFrame.from_dict(hits_df, orient='index', columns=fields)

    return hits_df

# RUNNING THE ANALYSIS

gbks = glob.glob(f"{GBK_DIR}/*.gbk")
proteins = extract_proteins_from_gbks(gbks)
print(len(proteins), "proteins extracted.")

hmms_path = f"hmms/*.hmm"  # HMM directory (from ktYPr repo)

proteins = extract_proteins_from_gbks(gbks)
results = run_variants(proteins, hmms_path)


# PRODUCE SUBTABLES for S17:

# --- Load the cutoff table
cutoffs_file = "ktypr/data/hmm_cutoffs_v20250704.tsv"
cutoffs = pd.read_csv(cutoffs_file, sep="\t")

# --- Assuming your results dataframe has a 'subject' column matching HMM_ID
# Add the cutoff values to results
results_with_cutoff = results.merge(cutoffs, left_on="subject", right_on="HMM_ID", how="left")

# --- Add a binary column: 1 if bitscore >= cutoff, else 0
results_with_cutoff["pass_cutoff"] = (results_with_cutoff["bitscore"] >= results_with_cutoff["Cutoff"]).astype(int)

# --- Optional: drop the HMM_ID column if redundant
results_with_cutoff = results_with_cutoff.drop(columns=["HMM_ID"])

# Add types
results_with_cutoff['type'] = [i.split('_')[0] for i in results_with_cutoff['query']]
results_with_cutoff['value'] = [int(i.split('_')[1]) for i in results_with_cutoff['query']]

results_with_cutoff.to_excel('./results.xlsx')

# And for the summary
import pandas as pd
import matplotlib.pyplot as plt

# Assume your dataframe is called df
summary = []

for locus, sub in results_with_cutoff.groupby("locus"):
    entry = {"locus": locus}
    
    # Mut: min value with pass_cutoff == 1
    mut_vals = sub[(sub["type"] == "Mut") & (sub["pass_cutoff"]==1)]["value"]
    entry["Mut_min"] = mut_vals.min() if not mut_vals.empty else None

    # N, C, Both: max value with pass_cutoff == 1
    for t in ["N", "C", "Both"]:
        vals = sub[(sub["type"] == t) & (sub["pass_cutoff"]==1)]["value"]
        entry[f"{t}_max"] = vals.max() if not vals.empty else None

    summary.append(entry)

summary_df = pd.DataFrame(summary)
summary_df.to_excel('./summary_results.xlsx', index=False)