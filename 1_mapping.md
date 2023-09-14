# Data mapping


Raw data processing was done with a pipeline I wrote here (https://github.com/xinsun1/tools_ngs/tree/main/tools_sunxin)

Here is the detail of our parameter choice
## 1. cutadapt
Double strand libary
```
cutadapt 
    -j #_threds
    -q 30,30
    -m 20
    -a N{3}AGATCGGAAGAGC
    -A N{3}AGATCGGAAGAGC
    -u 3
    -O 6
    -o trim_CMDARG4_1.fq -p trim_CMDARG4_2.fq
    ../0-raw/CMDARG5_1.fq.gz
    ../0-raw/CMDARG5_2.fq.gz
```

## 2. AdapterRemoval
We use AdapterRemoval to merge pairend read
```
AdapterRemoval
    --file1 trim_CMDARG1_1.fq
    --file2 trim_CMDARG1_2.fq
    --minalignmentlength 20
    --threads CMDARG2
    --collapse --basename CMDARG1
```

## 3. Mapping

```
bwa mem -t CMDARG1 CMDARG2 CMDARG3 CMDARG4
bwa aln -l 16500 -t CMDARG1 CMDARG2 CMDARG3
bwa samse CMDARG1 CMDARG2 CMDARG3
samtools view -b -h -q 20 CMDARG2 ' #[short_name, sam_f, type]
                            'CMDARG1_CMDARG3.sam -o CMDARG1_CMDARG3_f.bam
samtools sort -n -@ 20 -o CMDARG1_CMDARG2_f_nsort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f.bam
samtools fixmate -m CMDARG1_CMDARG2_f_nsort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f_nsort_fix.bam'
samtools sort -@ 20 -o CMDARG1_CMDARG2_f_sort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f_nsort_fix.bam
samtools stats CMDARG1_CMDARG2_f_sort.bam
samtools markdup -r CMDARG1_CMDARG2_f_sort.bam ' #[short_name, type]
                              'CMDARG1_CMDARG2_f_sort_rmdup.bam'
samtools index CMDARG1_CMDARG2_f_sort_rmdup.bam
```

Data from multiple libraries were then merged together and dedup again.
