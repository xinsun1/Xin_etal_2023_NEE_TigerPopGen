# GL replication of our results



## 1. error rate
```
#!/bin/bash


WDIR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/1-gl
REF=/share/Download_genomes/TigerRefGenome/P_tigris.scaffold.fa
BAM_LIST=$WDIR/list.bam

CHR_LIST=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/Cat_Chr_Map/autosome_chr_list


# gen anc
# cd $WDIR
# mkdir 0-angsd_fa
# cd 0-angsd_fa
# BAM=/share/users/sunx/modern_tiger_LYC/5.Tiger_aln_unique/RFET0002_clean_aln_pe_unique.bam
ID=RFET0002
#
# angsd -nThreads 4 -i $BAM -dofasta 2 -doCounts 1 -minQ 30 -minmapq 20 \
#       -setminDepthInd 3 -remove_bads 1 -uniqueOnly 1 \
#       -ref $REF -out $ID.minD3.q30.f2.cons
#
#

ANC=$WDIR/0-angsd_fa/$ID.minD3.q30.f2.cons.fa.gz


# err estimate
#### run per chr
#### set up by_chr_dir next time
cd $WDIR/3-qc
for i in $(less $CHR_LIST)
do
        RUN=0
        while (($RUN == 0));
        do
                PRUN=$(ps -ef | grep sunx | grep angsd | wc -l)

                if (($PRUN <= 20)); then
                        angsd -doAncError 1 -ref $ANC -anc $REF -out err_chr.$i -bam $BAM_LIST -C 50 -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -checkBamHeaders 0 -r ${i}: &
                        RUN=1
                else
                        sleep 30
                fi
        done
done
wait


#### merge result
cd $WDIR/3-qc
awk '{print "./by_chr/err_chr."$1".ancErrorChr"}' /share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/Cat_Chr_Map/autosome_chr_list > list.filename

awk 'FNR!=1 {for(i=1;i<=NF;i++){a[FNR-1][i]+=$i};tf=NF;tl=FNR-1;next}END{for(i=1;i<=tl;i++){for(j=1;j<=tf;j++){printf a[i][j]"\t";};print ""}}' $(< list.filename) > err_all.ancError.all


Rscript ../../program/angsd-master/R/estError.R file=err_all.ancError.all
```

## 2. GL call

```
#!/bin/bash

#### set working directory

WDIR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/1-gl
GL_FOLDER=$WDIR/gl_chr
REF=/share/Download_genomes/TigerRefGenome/P_tigris.scaffold.fa
BAM_LIST=$WDIR/list.bam
LIST_CHR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/Cat_Chr_Map/autosome_chr_list


cd $WDIR


#### GL call by chr

# if ! test -d gl_chr; then
#       mkdir gl_chr
# fi
#
# for i in $(less $LIST_CHR)
# do
#       RUN=0
#       while (($RUN == 0));
#       do
#               PRUN=$(ps -ef | grep sunx | grep angsd | wc -l)
#
#               if (($PRUN <= 20)); then
#                       angsd -GL 2 -out ./gl_chr/gl.$i -ref $REF -nThreads 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minMapQ 30 -minQ 20 -remove_bads 1 -uniqueOnly 1 -bam $BAM_LIST -r ${i}: &
#
#                       RUN=1
#               else
#                       sleep 30
#               fi
#       done
# done
#
# wait
#

#### merge
# exclude transversion


cd $WDIR
cd $GL_FOLDER

NSAMPLE=49

for i in $(less $LIST_CHR)
do
      # check maf file if -ref use $6,$8
      {
      zcat gl.${i}.mafs.gz | awk '{if(NR==1){print 1}else{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A") || $6 <0 || $6 >1 || $8/'$NSAMPLE' < 0 ){print 0}else{print 1} }}' > gl.${i}.is_tv_mafF_misF
      zcat gl.${i}.beagle.gz | awk 'NR==FNR {a[FNR]=$1;next}{if(a[FNR]==1){print $0}}' gl.${i}.is_tv_mafF_misF - > gl.${i}.tv_mafF_misF.beagle
      zcat gl.${i}.mafs.gz | awk 'NR==FNR {a[FNR]=$1;next}{if(a[FNR]==1){print $0}}' gl.${i}.is_tv_mafF_misF - > gl.${i}.tv_mafF_misF.m
afs
      } &
done
wait

CHR_H="scaffold1002"
head -1 gl.${CHR_H}.tv_mafF_misF.beagle >> $WDIR/gl_tv_mafF_misF.beagle
head -1 gl.${CHR_H}.tv_mafF_misF.mafs >> $WDIR/gl_tv_mafF_misF.mafs
for i in $(less $LIST_CHR)
do
       tail -n +2 gl.${i}.tv_mafF_misF.beagle >> $WDIR/gl_tv_mafF_misF.beagle
       tail -n +2 gl.${i}.tv_mafF_misF.mafs >> $WDIR/gl_tv_mafF_misF.mafs
done

gzip $WDIR/gl_tv_mafF_misF.beagle &
gzip $WDIR/gl_tv_mafF_misF.mafs &
wait
```

## 3. pca

```
#!/bin/bash

#### run
WDIR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/1-gl

cd $WDIR/1-pca

IN_BEAGLE=$WDIR/./gl_tv_mafF_misF.masked.beagle.gz
OUT="gl.m.tv.maf05"

for i in {2..10}
do
        pcangsd -b ${IN_BEAGLE} -e $i --iter 1000 -t 20 -o $OUT.e${i}.i1k
done
wait
```

Re-sampling

``` {r}
#### Downsize with procrustes projection ####
#### 1. random sampling per pop ----

#' generate random samples per pop
#' 
#' @param n number of samples per pop
#' @param df_sample df with three columns c(id(chr), pop(chr), for_sample(int, 1/0))
#' @param seed set seed for replication
#' 
#' @return two lists,  c(list_selected, list_projection)
#' 
#' @epport
#' 

generate_sampling <-  function(df_sample, n, seed=123){
    set.seed(seed)
    
    list_select = df_sample %>% 
        select(c(id,pop,for_sample)) %>% 
        filter(for_sample==1) %>% 
        group_by(pop) %>%
        slice_sample(n = n) %>% 
        pull(id)
    
    list_project = df_sample$id[! df_sample$id %in% list_select]
    
    return(list(list_select, list_project))
}

df_sampling = sample_meta %>% 
    mutate(for_sample=for_pca)

res_samping=generate_sampling(df_sampling, 3)
res_samping[[1]]
res_samping[[2]]

#### 2. do pca and projection ----

list_full = sample_meta %>%
    group_by(GL_order) %>%
    arrange(GL_order) %>%
    pull(id)
list_se = res_samping[[1]]
list_pro = res_samping[[2]]


pp_n3_pca12 =
    plot_pcangsd(
        './gl.tv.maf05.e6.i1k.cov', meta_legend, sample_meta,
        list_full = list_full, list_se = list_se,
        list_pro = list_pro, 
        list_id_plot = list_full,
        pro_dim = c(1,2)) +
    p_theme_pca +
    theme(
        # legend.position = c(0.8,0.7),
        legend.title = element_blank(),
        legend.background = element_rect(color="grey"),
        aspect.ratio = 1
    ) +
    guides(fill=guide_legend(ncol=1))

pp_pca12

#' generate pca by number of samples and seeds
#' 
#' @param list_n number of samples per pop
#' @param nrep number of replications per n
#' 
#' @return list_ggplots
#' 
#' @export
#' 

do_pca <- function(list_n, nrep, seed=12356){
    set.seed(seed)
    list_seed = runif(nrep * length(list_n)) * 1e6
    
    list_pca = list()
    idx_seed=1
    # run per n per rep
    for(i in 1:length(list_n)){
        for(j in 1:nrep){
            
            tmp_n = list_n[i]
            tmp_samping=generate_sampling(
                df_sampling, tmp_n,
                seed=list_seed[idx_seed])
            
            list_full = sample_meta %>%
                group_by(GL_order) %>%
                arrange(GL_order) %>%
                pull(id)
            list_se = tmp_samping[[1]]
            list_pro = tmp_samping[[2]]

            pp_n_pca12 =
                plot_pcangsd(
                    './gl.tv.maf05.e6.i1k.cov', meta_legend, sample_meta,
                    list_full = list_full, list_se = list_se,
                    list_pro = list_pro, 
                    list_id_plot = FALSE,
                    pro_dim = c(1,2)) +
                p_theme_pca +
                theme(
                    legend.position = "none",
                    legend.title = element_blank(),
                    legend.background = element_rect(color="grey"),
                    aspect.ratio = 1
                ) +
                labs(title =str_c('Samples per population:',tmp_n),
                     subtitle = str_c('Replication:',j))
            list_pca[[idx_seed]] = pp_n_pca12
            idx_seed = idx_seed + 1
        }
    }
    return(list_pca)
}

list_pca_plot = do_pca(c(3,4,5), 3)


# 3.1 merege multiple plots
library(cowplot)
# plot_grid(p_pca12, p_pca13, labels = NULL )
plot_grid(plotlist=list_pca_plot)

ggsave('pca_gl.grid.sample.rep.12x6.png',width = 12, height = 12,
       device = "png", dpi = 500)

```

## 4. ngsadm
```
#!/bin/bash


WDIR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/1-gl

# enter WDIR
cd $WDIR/2-ngsadm
IFS=$'\n'

IN_BGGZ=$WDIR/gl_tv_mafF_misF.masked.beagle.gz
BATCH='gl_m_tv_maf02_mis80'

# run rep
if ! test -d rep; then
    mkdir rep
fi
cd rep


for i in {2..9};
do

    for rep in {1..10}
    do

        RUN=0
        while (($RUN == 0));
        do
            PRUN=$(ps -ef | grep sunx | grep NGSadmix | wc -l)

            if (($PRUN <= 4)); then
                NGSadmix -likes $IN_BGGZ -outfiles ${BATCH}.${i}.${rep} -minMaf 0.02 -minInd 34 -P 10 -K $i &
                RUN=1
            else
                sleep 10
            fi
        done
    done
done
wait


# get best
cd $WDIR/2-ngsadm

if ! test -d best; then
        mkdir best
fi


for i in {2..9}
do
        NAME=$BATCH
        K=$i

        best_l=$(for rep in {1..10};do L=$(tail -1 ./rep/${NAME}.${K}.${rep}.log | awk '{split($2,a,"="); print a[2]}'); echo $rep $L; done | sort -k 2 -r -n | head -1 )
        best=$(echo $best_l | awk '{print $1}')
        best_like=$(echo $best_l | awk '{print $2}')

        cp ./rep/${NAME}.${K}.${best}.qopt ./best/${NAME}.${K}.Q
        echo $K $best $best_like >> ./best/best_summary.${NAME}

        # output detail
        for rep in {1..10};do L=$(tail -1 ./rep/${NAME}.${K}.${rep}.log | awk '{split($2,a,"="); print a[2]}'); echo $rep $L; done | sort -k 2 -r -n > ./best/best_detail.${NAME}.${K}
```


## 5. heterozygosity

```
#!/bin/bash

# set env

WDIR=/share/users/sunx/1-South_China_Tiger/2-holotype/0-AMO_sum/1-gl
LIST_NAME_BAM=$WDIR/list.id_bam.low
REF=/share/Download_genomes/TigerRefGenome/P_tigris.scaffold.fa

# enter WDIR
cd $WDIR/4-het
if ! test -d sfs_by_sample_tv;then
    mkdir sfs_by_sample_tv
fi
cd sfs_by_sample_tv


IFS=$'\n'

# check bam index
# check input bam file index
for i in $(less $LIST_NAME_BAM)
do
        ID=$(echo $i | awk '{print $1}')
        BAM=$(echo $i | awk '{print $2}')
        BAM_INDEX_1=$BAM".bai"
        BAM_INDEX_2=$(echo "$BAM" | sed -s 's/\.bam/\.bai/')

        if test -f $BAM_INDEX_1; then
                echo "bam index OK"
        else
                if test -f $BAM_INDEX_2; then
                        echo "bam index OK"
                else
                        module load samtools/1.10
                        samtools index $BAM &
                fi
        fi
done
wait


# get sfs
# switch to fatnode

for i in $(less $LIST_NAME_BAM);
do
    cd $WDIR/4-het/sfs_by_sample_tv

    ID=$(echo $i | awk '{print $1}')
        BAM=$(echo $i | awk '{print $2}')

    if ! test -d $ID; then
        mkdir $ID
    fi
    cd $ID

    RUN=0
    while (($RUN == 0));
    do
        PRUN=$(ps -ef | grep sunx | grep realSFS | wc -l)

        if (($PRUN <= 5)); then
            {
            /opt/install.bak/angsd.0.94/angsd/misc/realSFS -P 4 ${ID}.saf.idx > ${ID}.ml
            awk '{print $2/($1+$2)}' ${ID}.ml > ${ID}.het
            } &
            RUN=1
        else
            sleep 1
        fi
    done
    cd ..
done
wait
```
