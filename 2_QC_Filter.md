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





