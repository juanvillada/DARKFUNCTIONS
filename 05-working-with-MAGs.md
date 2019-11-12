# Working with MAGs

We will be working here with all MAGs (≥50% completeness, ≤10% contamination) from co-assemblies #1 and #2. We will start by renaming the contigs to make their names shorter and more informative adding the MAG they belong to.

### Get MAGs and rename contigs

```bash
cd $WORKDIR/05_MAGs/

conda activate anvio6

BINNING=01 # OR
BINNING=02

for CLUSTER in $(seq 51 100 | awk -v OFS='_' '{print "Bin", $0}'); do
  MAGs=$(ls $WORKDIR/04_BINNING/BINNING_"$BINNING"/CONCOCT_SPLIT/"$CLUSTER".SUMMARY/bin_by_bin | grep "Bin")

  for MAG in $MAGs; do
    anvi-script-reformat-fasta $WORKDIR/04_BINNING/BINNING_"$BINNING"/CONCOCT_SPLIT/"$CLUSTER".SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               -o REDUNDANT_MAGS/ASSEMBLY_"$BINNING"_$MAG.fa \
                               --simplify-names \
                               --prefix ASSEMBLY_"$BINNING"_$MAG
  done
done
```

---

Now we will run CHECKM and METAQUAST to compute some statistics about the MAGs.

### Analyse MAGs with CHECKM

```bash
cd $WORKDIR/05_MAGs

conda activate checkm

checkm tree REDUNDANT_MAGS \
            CHECKM \
            -x fa \
            -t 40

checkm tree_qa CHECKM \
               -f CHECKM/checkm_tree_qa.txt \
               -o 2 \
               --tab_table

checkm tree_qa CHECKM \
               -f CHECKM/checkm_tree.tre \
               -o 4

checkm lineage_set CHECKM \
                   CHECKM/lineage.ms

checkm analyze CHECKM/lineage.ms \
               c \
               CHECKM \
               -x fa \
               -t 40

checkm qa CHECKM/lineage.ms \
          CHECKM \
          -f CHECKM/checkm_qa.txt \
          -o 2 \
          -t 40 \
          --tab_table
```

### Plot CHECKM tree with ANVI'O

```bash
cd $WORKDIR/05_MAGs

conda activate anvio6

anvi-script-checkm-tree-to-interactive -t CHECKM/storage/tree/concatenated.tre \
                                       -o CHECKM_TREE

anvi-interactive -t CHECKM_TREE/newick.tree \
                 -d CHECKM_TREE/view_data.txt \
                 -p CHECKM_TREE/tree.db \
                 --manual
```

### Analyse MAGs with QUAST

```bash
cd $WORKDIR/05_MAGs

conda activate quast

quast.py REDUNDANT_MAGS/*.fa \
         -o QUAST \
         -t 20 \
         --no-plots
```

---

Here we will perform a prelimary phylogenomic analysis of the MAGs with GTDB-tk.

### Do phylogenomic analysis with GTDB-tk

```bash
cd $WORKDIR/05_MAGs

conda activate GTDBtk

gtdbtk classify_wf --genome_dir REDUNDANT_MAGS \
                   --out_dir GTDB \
                   --extension fa \
                   --cpus 40
```
