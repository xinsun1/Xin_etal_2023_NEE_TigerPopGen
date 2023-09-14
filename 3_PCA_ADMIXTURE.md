# PCA and ADMIXTURE

## 1. PCA
smartpca

```

```


## 2. ADMIXTURE
10 independent run for the following script

```
for i in {2..10}
do
    admixture --seed=$RANDOM --cv -j8 auto_neutral_noOUT_hq_trvn_recode12.ped ${i} > log_${i} 2>&1 &
done
```