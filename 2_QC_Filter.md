# QC Filter

Vcf file filter was performed using bcftools and custome scripts

`$TOOLS_SUNXIN` is https://github.com/xinsun1/tools_ngs/blob/main/tools_sunxin/


## 1. keep only bialleic SNPs
Keep only Scaffolds above 10KB, remove 
```
#!/bin/bash

bcftools filter -g 5 -T scaffold_min10kb -e 'BaseQRankSum<=-1.96 && BaseQRankSum>=-1.96 && MQRankSum <= -1.96 && F_MISSING > 0.5' 10th_raw.vcf  | bcftools view -m2 -M2 -v snps -O v -o 11th_hard.vcf
```

## 2. find neutral loci
```
#!/bin/bash
# filter for snps # ./filter_snp_al2.sh

vcf=$1

# bgzip -@ 80 ${vcf}
# bcftools index ${vcf}.gz

# generate auto, chrX, chrY
bcftools view -T ../../Cat_Chr_Map/x_chr_list.bed ${vcf} -O v -o chrX_${vcf} &
bcftools view -T ../../Cat_Chr_Map/y_chr_list.bed ${vcf} -O v -o chrY_${vcf} &
bcftools view -T ../../Cat_Chr_Map/autosome_chr_list.bed ${vcf} -O v -o auto_${vcf}

bcftools view -T ^../../Cat_Chr_Map/exon_1k_merge.bed chrX_${vcf} -O v -o chrX_neutral_${vcf} &
bcftools view -T ^../../Cat_Chr_Map/exon_1k_merge.bed chrY_${vcf} -O v -o chrY_neutral_${vcf} &
bcftools view -T ^../../Cat_Chr_Map/exon_1k_merge.bed auto_${vcf} -O v -o auto_neutral_${vcf}
```

## 3. filter for DP and MAC
Sample information for filter.
```
A3      Degraded        degraded        2.5
AMO3305 Degraded        degraded        3.5
AMO3308 Degraded        degraded        2.5
HNSD    Degraded        degraded        1.5
HPS     Degraded        degraded        5.5
NGH2    Degraded        degraded        7
PTV02   Degraded        degraded        5
PTV06   Degraded        degraded        8
PTV17   Degraded        degraded        4.5
RFET0002        Modern  fresh   12
RFET0003        Modern  fresh   12
RFET0005        Modern  fresh   12
RFET0006        Modern  fresh   12
RFET0007        Modern  fresh   12
RFET0011        Modern  fresh   12
RUSA21  Modern  degraded        8.5
SRR836361       OUT     fresh   26
SRR836372       OUT     fresh   26
Y14     Degraded        degraded        1.5
Y3      Degraded        degraded        1.5
Y7      Degraded        degraded        3
pti096  Modern  fresh   12
pti103  Modern  fresh   12
pti105  Modern  fresh   12
pti177  Modern  fresh   12
pti183  Modern  fresh   12
pti184  Modern  fresh   12
pti185  Modern  fresh   12
pti186  Modern  fresh   12
pti208  Modern  fresh   12
pti217  Modern  fresh   12
pti220  Modern  degraded        12
pti247  Modern  fresh   12
pti248  Modern  fresh   12
pti262  Modern  fresh   12
pti269  Modern  fresh   12
pti270  Modern  fresh   12
pti272  Modern  fresh   12
pti273  Modern  fresh   12
pti292  Modern  fresh   12
pti303  Modern  fresh   12
pti304  Modern  fresh   12
pti305  Modern  fresh   12
pti306  Modern  fresh   12
pti307  Modern  fresh   12
pti330  Modern  fresh   12
pti331  Modern  fresh   12
pti332  Modern  fresh   12
```

Filter was implemented with these two scripts
https://github.com/xinsun1/tools_ngs/blob/main/tools_sunxin/vcf_filter_new_step1.py
https://github.com/xinsun1/tools_ngs/blob/main/tools_sunxin/vcf_filter_new_step2.py

Thresholds are listed in the scripts.

## 4. SNP summary
SNP type summary from vcf files was conducted with this script
https://github.com/xinsun1/tools_ngs/blob/main/tools_sunxin/vcf_stat.py



