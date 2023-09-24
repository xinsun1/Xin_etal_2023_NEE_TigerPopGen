# Demographic modelling

## 1. PSMC

`run_psmc_pipe.sh`
```
#!/bin/bash

ID=$1

$TOOLS_SUNXIN/vcf2fa.py auto_neutral_withOUT_hqc.vcf ../P_tigris_mask_auto.fa ${ID} ${ID}_psmc 4

# generate PSMC fasta
/opt/install.bak/psmc-master/utils/fq2psmcfa -q 20 ${ID}_psmc.fa > ${ID}.psmcfa
/opt/install.bak/psmc-master/utils/splitfa ${ID}.psmcfa > ${ID}_split.psmcfa

# run PSMC
psmc -N25 -t20 -p "4+25*2+4+6" -o ${ID}.psmc ${ID}.psmcfa

# plot PSMC
psmc_plot.pl -g5 -u0.64e-8 plot_${ID} ${ID}.psmc
```
`run_psmc_bs.sh`
```
#!/bin/bash

bs_index=$1
sample=$2

psmc -N25 -t15 -b -p "4+25*2+4+6" -o ./psmc_bootstrap/${bs_index}.psmc ${sample}_split.psmcfa
```
For plot.
```
#!/bin/bash

for i in $(ls pti*psmc RFET*psmc)
do
    psmc_plot.pl -R -Y10 -g5 -u0.40e-8 plot_04_${i} ${i}
    awk '{print $1,$2,"m","'${i}'"}' OFS='\t' plot_04_${i}.0.txt >> ./psmc_plot_04_bs_txt
done

for i in $(ls PTV*psmc)
do
    psmc_plot.pl -R -Y10 -g5 -u0.40e-8 -M "PTV06=0.32" plot_04_${i} ${i}
    awk '{print $1,$2,"m","'${i}'"}' OFS='\t' plot_04_${i}.0.txt >> ./psmc_plot_04_bs_txt
done

for i in $(ls RUSA*psmc)
do
    psmc_plot.pl -R -Y10 -g5 -u0.40e-8 -M "RUSA21=0.2" plot_04_${i} ${i}
    awk '{print $1,$2,"m","'${i}'"}' OFS='\t' plot_04_${i}.0.txt >> ./psmc_plot_04_bs_txt
done

for i in $(ls HPS*psmc)
do
    psmc_plot.pl -R -Y10 -g5 -u0.40e-8 -M "HPS=0.32" plot_04_${i} ${i}
    awk '{print $1,$2,"m","'${i}'"}' OFS='\t' plot_04_${i}.0.txt >> ./psmc_plot_04_bs_txt
done
```
Plot with R.
``` {r}
setwd('~/Documents/Projects/Tiger_ancient_historical/South_China_tiger/Progress_project/Demographic_history/191125_psmc/')


# read raw file
raw = read.table('psmc_plot_bs_txt', header = F, stringsAsFactors = F)
raw = read.table('psmc_plot_04_bs_txt', header = F, stringsAsFactors = F)


raw_new = c()
for(i in seq(1,nrow(raw))){
    if(raw[i,1]!=0){
        tmp_row = raw[i-1,]
        tmp_row[,1] = raw[i,1]
        raw_new = rbind(raw_new, tmp_row)
        raw_new = rbind(raw_new, raw[i,])
    }
    else{
        raw_new = rbind(raw_new, raw[i,])
    }
}

write.table(raw_new, file='psmc_plot_bs_txt_len', quote = FALSE, sep = '\t',
            col.names = FALSE, row.names = FALSE)


library(ggplot2)

colnames(raw) = c('time', 'Ne', 'type', 'indi')


p_theme = theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
                panel.grid.major = element_line(color = NA),
                panel.grid.minor = element_line(color = NA),
                #axis.text.x = element_text(angle = 90),
                axis.title.y = element_blank(), axis.title.x = element_blank(),
                legend.key = element_rect(colour = NA, fill = NA),
                legend.position = c(0.3,0.7),
                legend.background = element_rect(fill=NA)
                
)


raw_m = raw[raw$type=='m',]
type_raw = unlist(lapply(raw_m$indi,function(x)(unlist(strsplit(x, '\\.'))[[1]])))
raw_m = cbind(raw_m, type_raw)
colnames(raw_m) = c('time', 'Ne', 'type', 'indi', 'type_indi')
raw_m = raw_m[!type_raw %in% c('pti217', 'pti220'),]

raw_bs = raw[raw$type=='bs',]
raw_bs = cbind(raw_bs, raw_bs$type)
colnames(raw_bs) = c('time', 'Ne', 'type', 'indi', 'type_indi')
type_indi_bs = unlist(lapply(raw_bs$indi,function(x)(unlist(strsplit(x, '_'))[[1]])))
raw_bs = raw_bs[!type_indi_bs %in% c('pti217', 'pti220'),]

raw_order = rbind(raw_bs, raw_m)

color_panel = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',
                '#ffff33','#a65628','#f781bf','#999999')
raw_order_10k = raw_order[raw_order$time>=1e3,]

ggplot(raw_order_10k) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + p_theme +
    geom_step(data = subset(raw_order_10k, type == 'bs'),
              aes(x=time, y=Ne, group=indi), color='grey', alpha=0.5) +
    geom_step(data = subset(raw_order_10k, type == 'm'),
              aes(x=time, y=Ne, group=indi, color=type_indi), size= 1) +
    scale_color_manual(values = color_panel) + annotation_logticks(sides = 'tb') +
    ylim(0,20)


# test plot 
ggplot(raw_order) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + p_theme +
    geom_step(data = subset(raw_order, type == 'm'),
              aes(x=time, y=Ne, group=indi, color=type_indi), size= 1) +
    scale_color_manual(values = color_panel) + annotation_logticks(sides = 'tb') +
    ylim(0,20)


# read climate data
clim = read.table('climate_1999_no_head', header = FALSE, stringsAsFactors = FALSE)
clim = clim[,c(2,4)]
colnames(clim) = c('years', 'temperature')

ggplot(clim[clim$years >= 1e4,]) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + p_theme + geom_line(aes(x=years, y=temperature)) +
     annotation_logticks(sides = 'tb') +
    scale_y_continuous(breaks = c(seq(-10,4,2)))

clim_mod = cbind(clim, rep('t', nrow(clim)),
                 rep('t', nrow(clim)),
                 rep('t', nrow(clim)),
                 rep('t', nrow(clim)))
colnames(clim_mod) = c('time', 'Ne', 'type', 'indi', 'type_indi', 'facet')
clim_mod = clim_mod[clim_mod$time>=1e4,]

raw_facet = cbind(raw_order_10k, rep('psmc', nrow(raw_order_10k)))
colnames(raw_facet) = c('time', 'Ne', 'type', 'indi', 'type_indi', 'facet')
psmc_clim = rbind(raw_facet, clim_mod)


ggplot(psmc_clim) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + p_theme +
    geom_step(data = subset(psmc_clim, type == 'bs'),
              aes(x=time, y=Ne, group=indi), color='grey', alpha='0.5') +
    geom_step(data = subset(psmc_clim, type == 'm'),
              aes(x=time, y=Ne, group=indi, color=type_indi), size= 1) +
    geom_line(data = subset(psmc_clim, type == 't'),
              aes(x=time, y=Ne)) +
    facet_grid(facet ~ ., scales = "free_y") + 
    scale_color_manual(values = color_panel) + annotation_logticks(sides = 'tb')
```

## 2. BEAST

xml file `./files/50M_all_nop_cep_node`

run replications
```
#!/bin/bash

XML=$1

for i in {1..4}
do
    mkdir ${i}-rep${i}
    cp $XML ./${i}-rep${i}
    cd ./${i}-rep${i}
    beast2 $XML > run_log_${XML}.log 2>&1 &
    cd ..
done
```

## 3. GPhoCS

Generate loci
```
#!/bin/bash

$TOOLS_SUNXIN/gen_loci_window.py -f name_for_gen_gphocs -w 1000 -s 51000 -l 10 --chr $1 -n 0.01 -o gphocs_loci
```

Control files in `./files/control_gphocs`

Pilot run on randomly selected 5k loci with 3 replication.
Final run with all loci with 2 replications.

Process tracer results
``` {r}
# '''
# read tracer statistics 
# '''

wdir = '/Users/sunxin/Documents/ancient_tiger/South_China_tiger/Progress_project/Population_history_modeling/190924_gphocs_pilot/'
setwd(wdir)

# read tracer file 

read_trace = function(file_name, raw, tmp_var_list){
    fh = read.table(file_name, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
    rownames(fh) = fh[,1]
    fh = fh[,2:ncol(fh)]
    
    for(i in 1:ncol(fh)){
        # generate variable information
        col_list = unlist(strsplit(fh[1,i], ': '))
        var_name = col_list[2]
        var_name_list = unlist(strsplit(var_name, '_'))
        var_type = var_name_list[1]
        var_info = paste(var_name_list[2:length(var_name_list)], collapse="_")
        
        run_name_all = col_list[1]
        run_name_list = unlist(strsplit(run_name_all, '_'))
        run_name_set = run_name_list[2]
        run_name_info = paste(run_name_list[3:length(run_name_list)], collapse="_")
        run_name_batch = paste(c(run_name_set, run_name_info, var_name), collapse="_")
        if(run_name_batch %in% tmp_var_list){
            run_name_rep = 'rep2'
        }else{
            run_name_rep = 'rep1'
            tmp_var_list = c(tmp_var_list, run_name_batch)
        }
        run_batch_mark = paste(c(run_name_rep, run_name_set, run_name_info), collapse="_")
        # add mean value
        mean_value = as.numeric(fh[2,i])
        # add 95% HPD
        hpd_up = as.numeric(strsplit(strsplit(fh[9,i], ", ")[[1]][2], '[]]')[[1]][1])
        hpd_down = as.numeric(strsplit(strsplit(fh[9,i], ", ")[[1]][1], '[[]')[[1]][2])
        raw = rbind(raw,
                    c(mean_value, hpd_down, hpd_up,
                      run_name_rep, run_name_set, run_name_info,
                      var_name, var_type, var_info,
                      run_batch_mark))
    }
    out_list = list(raw, tmp_var_list)
    return(out_list)
}
    # test read_trace
out = read_trace('result_20w_mig', raw=raw, tmp_var_list = tmp_var_list)
raw = out[[1]]
tmp_var_list = out[[2]]


# read in all trace file
raw = c()   # dataframe for store
tmp_var_list = c()
file_list = list.files()
file_list = file_list[grep('result', file_list)]
file_list = c('Result_20w_mig', 'Result_20w_theta_tau')



for(i in 1:length(file_list)){
    out = read_trace(file_list[i], raw=raw, tmp_var_list = tmp_var_list)
    raw = out[[1]]
    tmp_var_list = out[[2]]
}
colnames(raw) = c('mean', 'hpd_down', 'hpd_up',
                  'loci_set', 'indi_set', 'mig_set',
                  'var', 'var_type', 'var_info',
                  'run_id')

# plot theta and tau
library(ggplot2)

data_frame_raw = as.data.frame(raw, stringsAsFactors = FALSE)
# combine to NA
data_frame_raw$indi_set[is.na(data_frame_raw$indi_set)] = 'combine'
theta_raw = data_frame_raw[data_frame_raw$var_type=='theta',]
tau_raw = data_frame_raw[data_frame_raw$var_type=='tau',]



    # run_rep loci set
    # run_data_set individual set
p_theme = theme(panel.background = element_rect(fill = NA, colour = "black", linetype = 1),
                panel.grid.major = element_line(color = NA),
                panel.grid.minor = element_line(color = NA),
                #axis.text.x = element_text(angle = 90),
                axis.title.y = element_blank(), axis.title.x = element_blank(),
                legend.position = "none",
                strip.background=element_blank()
)

color_panel = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',
                '#ffff33','#a65628','#f781bf','#999999')
    # theta plot
ggplot(data = theta_raw) + geom_point(aes(x=as.factor(var),
                                        y=as.numeric(mean),
                                        group=as.factor(run_id),
                                        color=as.factor(indi_set),
                                        shape=as.factor(loci_set)),
                                      position=position_dodge(width = 0.90)) +
    geom_errorbar(aes(x=as.factor(var),
                      ymin=as.numeric(hpd_down),
                      ymax=as.numeric(hpd_up),
                      group=as.factor(run_id),
                      color=as.factor(indi_set)),
                  position=position_dodge(width = 0.90))
    # theta plot facet
ggplot(data = theta_raw) + geom_point(aes(x=as.numeric(as.factor(run_id)),
                                          y=as.numeric(mean),
                                          color=as.factor(indi_set),
                                          shape=as.factor(loci_set)),
                                      position=position_dodge(width = 0.90)) +
    geom_errorbar(aes(x=as.numeric(as.factor(run_id)),
                      ymin=as.numeric(hpd_down),
                      ymax=as.numeric(hpd_up),
                      color=as.factor(indi_set)),
                  position=position_dodge(width = 0.90)) +
    facet_wrap(.~as.factor(var), scales="free_y") +
    p_theme + scale_color_manual(values=color_panel)


    # tau plot facet
ggplot(data = tau_raw) + geom_point(aes(x=as.numeric(as.factor(run_id)),
                                          y=as.numeric(mean),
                                          color=as.factor(indi_set),
                                          shape=as.factor(loci_set)),
                                      position=position_dodge(width = 0.90)) +
    geom_errorbar(aes(x=as.numeric(as.factor(run_id)),
                      ymin=as.numeric(hpd_down),
                      ymax=as.numeric(hpd_up),
                      color=as.factor(indi_set)),
                  position=position_dodge(width = 0.90)) +
    facet_wrap(.~as.factor(var), scales="free_y") +
    p_theme + scale_color_manual(values=color_panel)

# generate total mig rate
mig_raw = data_frame_raw[data_frame_raw$var_type=='m',]


mig_ref_mig = c(c('SUM','TIG'), c('SUM','RUSA'), c('SUM','AMO'), c('SUM','PTV'),
                c('SUM','AMUR'), c('SUM','COR'), c('SUM','JAX'),
                c('TIG','RUSA'), c('TIG','AMO'), c('TIG','PTV'), c('TIG','AMUR'),
                c('TIG','COR'), c('TIG','JAX'),
                c('RUSA','AMO'), c('RUSA','PTV'), c('RUSA','AMUR'), c('RUSA','COR'),
                c('RUSA','JAX'),
                c('AMO','PTV'), c('AMO','AMUR'), c('AMO','COR'), c('AMO','JAX'),
                c('PTV','AMUR'), c('PTV','COR'), c('PTV','JAX'), 
                c('AMUR','COR'), c('AMUR','JAX'), 
                c('COR','JAX'),
                c('AN_AP', 'AN_CJ'), c('AN_AP', 'RUSA')
)
mig_ref_mig = matrix(mig_ref_mig, nrow=2)
mig_ref_tau = c('tau_AN_TIG_div', 'tau_AN_RA', 'tau_AN_AP', 'tau_AN_AMP',
  'tau_AN_AMP', 'tau_AN_CJ', 'tau_AN_CJ',
  'tau_AN_RA', 'tau_AN_AP', 'tau_AN_AMP', 'tau_AN_AMP',
  'tau_AN_CJ', 'tau_AN_CJ',
  'tau_AN_AP', 'tau_AN_AMP', 'tau_AN_AMP', 'tau_AN_CJ',
  'tau_AN_CJ',
  'tau_AN_AMP', 'tau_AN_AMP', 'tau_AN_CJ', 'tau_AN_CJ',
  'tau_AN_AMP', 'tau_AN_CJ', 'tau_AN_CJ',
  'tau_AN_CJ', 'tau_AN_CJ',
  'tau_AN_CJ',
  'tau_AN_RA,tau_AN_AP', 'tau_AN_RA,tau_AN_AP'
)

mig_ref = cbind(mig_ref_mig[1,],
                mig_ref_mig[2,],
                mig_ref_tau)
mig_ref_band = apply(mig_ref, 1, function(x){return(c(paste(x[c(1,2)], collapse="->"),
                                                      paste(x[c(2,1)], collapse="->")))})
mig_ref = cbind(mig_ref, t(mig_ref_band))

total_mig_rate = c()
for(i in 1:nrow(mig_raw)){
    tmp_mig_info = mig_raw$var_info[i]
    tmp_mig_runid = mig_raw$run_id[i]
    if(tmp_mig_info %in% mig_ref[,4]){
        tau_tmp = mig_ref[mig_ref[,4]==tmp_mig_info,3]
    }else{
        tau_tmp = mig_ref[mig_ref[,5]==tmp_mig_info,3]
    }
    tau_tmp_list = unlist(strsplit(tau_tmp, ","))
    if(length(tau_tmp_list)==1){
        tau_tmp_mig = tau_raw$mean[tau_raw$run_id==tmp_mig_runid & tau_raw$var==tau_tmp_list[1]]
        tau_tmp_mig = as.numeric(tau_tmp_mig)
    }else{
        tau_tmp_mig_1 = tau_raw$mean[tau_raw$run_id==tmp_mig_runid & tau_raw$var==tau_tmp_list[1]]
        tau_tmp_mig_2 = tau_raw$mean[tau_raw$run_id==tmp_mig_runid & tau_raw$var==tau_tmp_list[2]]
        tau_tmp_mig = abs(as.numeric(tau_tmp_mig_1) - as.numeric(tau_tmp_mig_2))
    }
    # theta_tau_print 10000
    # mig_print     0.001
    tmp_total_mig = as.numeric(mig_raw$mean[i]) * tau_tmp_mig / 10000 / 0.001
    tmp_total_mig_max = as.numeric(mig_raw$hpd_up[i]) * tau_tmp_mig / 10000 / 0.001
    tmp_total_mig_min = as.numeric(mig_raw$hpd_down[i]) * tau_tmp_mig / 10000 / 0.001
    
    total_mig_rate = rbind(total_mig_rate,
                           c(tmp_total_mig, tmp_total_mig_min, tmp_total_mig_max))
}
colnames(total_mig_rate) = c('total_mig_rate', 'total_mig_rate_min', 'total_mig_rate_max')

mig_raw_with_rate = cbind(mig_raw, total_mig_rate) 

ggplot(data = mig_raw_with_rate) + geom_point(aes(x=as.numeric(as.factor(run_id)),
                                        y=as.numeric(total_mig_rate),
                                        color=as.factor(indi_set),
                                        shape=as.factor(loci_set)),
                                    position=position_dodge(width = 0.90)) +
    geom_errorbar(aes(x=as.numeric(as.factor(run_id)),
                      ymin=as.numeric(total_mig_rate_min),
                      ymax=as.numeric(total_mig_rate_max),
                      color=as.factor(indi_set)),
                  position=position_dodge(width = 0.90)) +
    facet_wrap(.~as.factor(var), scales="free_y") +
    p_theme + scale_color_manual(values=color_panel)


ggplot(data = mig_raw_with_rate) + geom_point(aes(x=as.numeric(as.factor(run_id)),
                                                  y=as.numeric(total_mig_rate),
                                                  color=as.factor(indi_set),
                                                  shape=as.factor(loci_set)),
                                              position=position_dodge(width = 0.90)) +
    geom_errorbar(aes(x=as.numeric(as.factor(run_id)),
                      ymin=as.numeric(total_mig_rate_min),
                      ymax=as.numeric(total_mig_rate_max),
                      color=as.factor(indi_set)),
                  position=position_dodge(width = 0.90)) +
    facet_wrap(.~var, scales="free_y") +
    geom_hline(yintercept=0.01, linetype=2) +
    p_theme + scale_color_manual(values=color_panel)
    


#### output table ####
write.table(tau_raw,file = 'tau_all', quote = F, sep = '\t')
write.table(theta_raw,file = 'theta_all.txt', quote = F, sep = '\t')
write.table(mig_raw_with_rate,file = 'mig_all.txt', quote = F, sep = '\t')

```
