# Extract nt and aa in two files
import glob
from Bio import SeqIO
import os

def extract_nucleotides_from_gbks(gbk_files):
    nucleotides = {}

    for gbk in gbk_files:
        kID = os.path.basename(gbk).split('/')[-1].split('_')[0]
        for record in SeqIO.parse(gbk, "genbank"):
            for feat in record.features:
                if feat.type == "CDS":
                    locus = feat.qualifiers.get("gene", [None])[0]

                    if locus and locus not in nucleotides:
                        seq = feat.location.extract(record.seq)
                        nucleotides[kID+'__'+locus] = str(seq)

    return nucleotides

def extract_proteins_from_gbks(gbk_files):
    proteins = {}
    for gbk in gbk_files:
        kID = os.path.basename(gbk).split('/')[-1].split('_')[0]
        for record in SeqIO.parse(gbk, "genbank"):
            for feat in record.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    locus = feat.qualifiers.get("gene", [None])[0]
                    seq = feat.qualifiers["translation"][0]
                    if locus and locus not in proteins:
                        proteins[kID+'__'+locus] = seq
    return proteins

def write_fasta(seqs, out):
    with open(out, "w") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")

# Gbk files path
gbk_files = glob.glob("ktypr/data/reference_clusters/*.gbk")

aa = extract_proteins_from_gbks(gbk_files)
nt = extract_nucleotides_from_gbks(gbk_files)

write_fasta(aa, "./blast/proteins.faa")
write_fasta(nt, "./blast/nucleotides.fna")