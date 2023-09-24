# PCA and ADMIXTURE

## 1. PCA
smartpca

```
genotypename:    auto_neutral_noOUT_hq_trvn.geno
snpname:         auto_neutral_noOUT_hq_trvn.snp
indivname:       auto_neutral_noOUT_hq_trvn_pop.indi
evecoutname:     auto_neutral_noOUT_hq_trvn.evec
evaloutname:     auto_neutral_noOUT_hq_trvn.eval
outliermode:     2
fastmode:        NO
numthreads:      40
numchrom:        1
poplistname:     pop_list
```

## 2. ADMIXTURE
10 independent run for the following script

```
for i in {2..10}
do
    admixture --seed=$RANDOM --cv -j8 auto_neutral_noOUT_hq_trvn_recode12.ped ${i} > log_${i} 2>&1 &
done
```


## 3. Phylogeny
RAxML 
run with default setting.
```
raxmlHPC-HYBRID-AVX -p 21370923 -x 1123840 -T 40 -# 100 -m GTRGAMMA -s ../../auto_neutral_withOUT_psdo_trvn.phylip -n auto_neutral_withOUT_psdo_trvn
```

Treemix
```
#!/bin/bash

treemix=$1

n=1
for i in {0..9}
do
    mkdir rep_${i}
    cd rep_${i}
    for j in {1..100}
    do
        if ((i==0)); then
            treemix -seed $RANDOM -i ../${treemix} -root PLE,PUN -k 1000 -o out_k1k_m${i}_rep${j} > log_k1k_${i}_${j} &
        else
            treemix -seed $RANDOM -i ../${treemix} -root PLE,PUN -k 1000 -m ${i} -o out_k1k_m${i}_rep${j} > log_k1k_${i}_${j} &
        fi
        ((n+=1))
        echo $n
        if ((n==40)); then
            sleep 600
            n=0
        fi
    done
    cd ..
done
```