# Comparative with Kaptive 

- Downloaded the database from https://github.com/rgladstone/EC-K-typing?tab=readme-ov-file
- Installed kaptive from https://kaptive.readthedocs.io/en/latest/Installation.html

Then run:

```bash
conda activate kaptive

kaptive assembly ./EC-K-typing-main/DB/EC-K-typing_group2and3_v3.0.0.gbk ./kref_genomes/*.fasta -o ./comp_data/kaptive_kref.tsv -j ./comp_data/jsons/kaptive_kref.json

kaptive assembly ./EC-K-typing-main/DB/EC-K-typing_group2and3_v3.0.0.gbk ./kref10/01_fasta/*.fna -o ./comp_data/kaptive_kref10.tsv -j ./comp_data/jsons/kaptive_kref10.json

kaptive assembly ./EC-K-typing-main/DB/EC-K-typing_group2and3_v3.0.0.gbk ./ktypr_reference/*.fasta -o ./comp_data/kaptive_bioinforef.tsv -j ./comp_data/jsons/kaptive_bioinforef.json

kaptive assembly ./EC-K-typing-main/DB/EC-K-typing_group2and3_v3.0.0.gbk ./EC-K-typing-main/DB/split_fasta/*.fasta -o ./comp_data/kaptive_kaptive.tsv -j ./comp_data/jsons/kaptive_kaptive.json
```

Kaptive with kTYPr custom DB (download (here)[https://github.com/SushiLab/kTYPr/blob/master/kaptive_db/KapsDB_for_Kaptive.gbk]):

```bash
conda activate kaptive

kaptive assembly ./KapsDB_for_Kaptive.gbk ./kref/kref_genomes/*.fasta -o ./comp_data/kaptive_ktyprDB_kref.tsv -j ./comp_data/jsons/kaptive_ktyprDB_kref.json

kaptive assembly ./KapsDB_for_Kaptive.gbk ./kref10/01_fasta/*.fna -o ./comp_data/kaptive_ktyprDB_kref10.tsv -j ./comp_data/jsons/kaptive_ktyprDB_kref10.json

kaptive assembly ./KapsDB_for_Kaptive.gbk ./ktypr_reference/*.fasta -o ./comp_data/kaptive_ktyprDB_bioinforef.tsv -j ./comp_data/jsons/kaptive_ktyprDB_bioinforef.json

kaptive assembly ./KapsDB_for_Kaptive.gbk ./EC-K-typing-main/DB/split_fasta/*.fasta -o ./comp_data/kaptive_ktyprDB_kaptive.tsv -j ./comp_data/jsons/kaptive_ktyprDB_kaptive.json
```

For ktypr on kaptive

```bash
ktypr -i ./EC-K-typing-main/DB/split_gbk/ -o ./comp_data/kaptive_ktypr/ --mode 1 -n -1 -v -p ktypr_wg -c
ktypr -i./ktypr_reference/ -o ./comp_data/ktypr_ktypr/ --mode 1 -n -1 -v -c
```