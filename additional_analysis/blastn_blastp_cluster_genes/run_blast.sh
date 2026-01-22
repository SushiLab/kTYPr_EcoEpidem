cd ./blast/

makeblastdb -in nucleotides.fna -dbtype nucl -out nt_db
makeblastdb -in proteins.faa    -dbtype prot -out aa_db

HEADER="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

blastn -query nucleotides.fna -db nt_db -out tmp_nt.tsv \
       -outfmt 6 -num_threads 8

echo -e "$HEADER" | cat - tmp_nt.tsv > nt_vs_nt.tsv
rm tmp_nt.tsv

blastp -query proteins.faa -db aa_db -out tmp_aa.tsv \
       -outfmt 6 -num_threads 8

echo -e "$HEADER" | cat - tmp_aa.tsv > aa_vs_aa.tsv
rm tmp_aa.tsv