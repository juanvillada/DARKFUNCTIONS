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

checkm lineage_wf REDUNDANT_MAGS \
                  CHECKM \
                  -x fa \
                  -t $SLURM_CPUS_PER_TASK
```
