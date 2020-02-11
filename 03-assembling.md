# Metagenome assembling

Samples appear too different for a good co-assembly. We will then do two co-assemblies, one with the bog samples only and the other with the remaining samples. This is done by assigning the variables $ASSEMBLY, $R1 and $R2. We will use MEGAHIT as it is fast and has been shown to produce good assemblies (e.g. Quince et al. 2017, Sczyrba et al. 2017, van der Walt et al. 2017).

### Define assembly and create list of file names

```bash
ASSEMBLY=ASSEMBLY_01
R1=$(cat $WORKDIR/SAMPLES.txt | egrep -v [mo]122 | awk -v ORS="," '{print "../01_TRIMMED_DATA/" $0 "_R1_trimmed.fastq"}' | sed 's/,$/\n/')
R2=$(cat $WORKDIR/SAMPLES.txt | egrep -v [mo]122 | awk -v ORS="," '{print "../01_TRIMMED_DATA/" $0 "_R2_trimmed.fastq"}' | sed 's/,$/\n/')

#OR

ASSEMBLY=ASSEMBLY_02
R1=$(cat $WORKDIR/SAMPLES.txt | egrep [mo]122 | awk -v ORS="," '{print "../01_TRIMMED_DATA/" $0 "_R1_trimmed.fastq"}' | sed 's/,$/\n/')
R2=$(cat $WORKDIR/SAMPLES.txt | egrep [mo]122 | awk -v ORS="," '{print "../01_TRIMMED_DATA/" $0 "_R2_trimmed.fastq"}' | sed 's/,$/\n/')
```

### Assemble reads with MEGAHIT

```bash
cd $WORKDIR/03_ASSEMBLIES

module load intelmpi
module load mkl
module load python/3.5.3
module load megahit

megahit -1 $(echo $R1) \
        -2 $(echo $R2) \
        --out-dir $ASSEMBLY \
        --min-contig-len 1000 \
        --k-min 57 \
        --k-max 157 \
        --k-step 10 \
        --memory 0.8 \
        --num-cpu-threads 40
```

### Check assemblies with METAQUAST

```bash
cd $WORKDIR/03_ASSEMBLIES

module load biokit

metaquast.py ASSEMBLY_*/final.contigs.fa \  
             -o ASSEMBLIES_QC \  
             -t 16 \  
             --fast
```
