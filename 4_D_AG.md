# D-statistics and AdmixtureGraph

## 1. D-statistics
D-statistics was run with admixtools2.
``` {r}
# Final Dstat for tiger manuscript

# 1. set env --------------------------------------------------------------
library(tidyverse)
library(admixtools)
setwd("/Users/sunxin/Documents/Projects/2019_Tiger_ancient_historical/South_China_tiger/Progress_project/3.Gene_flow/2023_adm2/")

# 2. load run_matrix and run ----------------------------------------------
## 2.1 run previous Dstat ----
run_pops = read_delim(
    file = './D_pop_indi_list_run',
    delim = '\t',
    col_names = c('pop1', 'pop2', 'pop3', 'pop4'))

res_d = qpdstat("auto_neutral_withOUT_psdo_trvn", pop1 = run_pops,
                f4mode = FALSE) 
save(res_d, file = "res_d.rdata")

## 2.2 run added Dstat ----
# D(OUT, SUM, TIG, X)
# D(COR/JAX, AMO_1, AMO_X)
# D(O,TIG, COR/JAX, XX)

list_popX = sample_meta %>%
    dplyr::select(id) %>%
    mutate(id=case_when(
        id=="Ptip121" ~ "PtiP121",
        id=="Pti217" ~ "pti217",
        TRUE ~ id)) %>% 
    dplyr::pull()

prefix="auto_neutral_withOUT_psdo_trvn"

test = tibble(
    pop1=list(c("SRR836361", "SRR836372")),
    pop2=list(c("pti305","pti269")),
    pop3=list("HPS"),
    pop4=list(list_popX))

run_matrix = bind_rows(
    # D(OUT, SUM, TIG, X)
    tibble(
        pop1=list("SRR836361", "SRR836372"),
        pop2=list("pti096"),
        pop3=list("pti103"),
        pop4=list(list_popX)) %>%
        unnest() %>%
        expand(pop1, pop2, pop3, pop4),
    # D(OUT, COR/JAX, AMO_1, AMO_X)
    tibble(
        pop1=list("SRR836361", "SRR836372"),
        pop2=list("pti305","pti269"),
        pop3=list("HPS"),
        pop4=list(list_popX)) %>% 
        unnest() %>%
        expand(pop1, pop2, pop3, pop4),
    # D(O,TIG, COR/JAX, X)
    tibble(
        pop1=list("SRR836361", "SRR836372"),
        pop2=list("pti103"),
        pop3=list("pti305","pti269"),
        pop4=list(list_popX)) %>%
        unnest() %>%
        expand(pop1, pop2, pop3, pop4),
)

res_d = qpdstat(
    prefix,
    pop1 = run_matrix,
    f4mode = FALSE
)  

save(res_d, file="res.D.add_b1.rdata")

## 2.3 run f3 ----
# f3(out, RUSA/VIR, X)
list_popX = sample_meta %>%
    dplyr::select(id) %>%
    mutate(id=case_when(
        id=="Ptip121" ~ "PtiP121",
        id=="Pti217" ~ "pti217",
        TRUE ~ id)) %>% 
    dplyr::pull()

prefix="auto_neutral_withOUT_psdo_trvn"

pop1=c("SRR836372", "SRR836361")
pop2=c("RUSA21","PTV02", "PTV06", "PTV17")
pop3=c(list_popX)

res_f3 = qp3pop(
    prefix,
    pop1 = pop1,
    pop2 = pop2,
    pop3 = pop3,
    auto_only=FALSE,
    outgroupmode=TRUE,
    apply_corr=FALSE, # test with TRUE or FALSE, small sample size correction
    verbose=TRUE
)
save(res_f3, file="res.f3.rusa_vir.rdata")

## 2.2 run added Dstat2 ----
# D(X, AMO_1, AMO_X)
# D()
list_popX = sample_meta %>%
    dplyr::select(id) %>%
    mutate(id=case_when(
        id=="Ptip121" ~ "PtiP121",
        id=="Pti217" ~ "pti217",
        TRUE ~ id)) %>% 
    dplyr::pull()

prefix="auto_neutral_withOUT_psdo_trvn"

run_matrix = bind_rows(
    # D(O,TIG, COR/JAX, X)
    tibble(
        pop1=list(c("SRR836361", "SRR836372")),
        pop2=list("pti103"),
        pop3=list(c("pti305","pti269")),
        pop4=list(list_popX)) %>%
        unnest(pop1) %>% unnest(pop2) %>% 
        unnest(pop3) %>% unnest(pop4),
    # D(OUT, X, AMO_1, AMO_X)
    tibble(
        pop1=list(c("SRR836361", "SRR836372")),
        pop2=list(c("pti305","pti269", "pti096", "pti103", "RFET0003", "PTV06")),
        pop3=list("HPS"),
        pop4=list(list_popX)) %>% 
        unnest(pop1) %>% unnest(pop2) %>% 
        unnest(pop3) %>% unnest(pop4),
)

res_d = qpdstat(
    prefix,
    pop1 = run_matrix,
    f4mode = FALSE
)  

save(res_d, file="res.D.add_b2.rdata")

## 2.2 run added Dstat3 ----
# D(X, AMO/ALT/VIR, RUSA)
# D()
list_popX = sample_meta %>%
    dplyr::select(id) %>%
    mutate(id=case_when(
        id=="Ptip121" ~ "PtiP121",
        id=="Pti217" ~ "pti217",
        TRUE ~ id)) %>% 
    dplyr::pull()

prefix="auto_neutral_withOUT_psdo_trvn"

run_matrix = bind_rows(
    # D(X, AMO/ALT/VIR, RUSA)
    tibble(
        pop1=list(c("SRR836361", "SRR836372")),
        pop2=list(list_popX),
        pop3=list("HPS", "RFET0003", "PTV06"),
        pop4=list("RUSA21") ) %>% 
        unnest(pop1) %>% unnest(pop2) %>% 
        unnest(pop3) %>% unnest(pop4),
    
    
)

res_d = qpdstat(
    prefix,
    pop1 = run_matrix,
    f4mode = FALSE
)  

save(res_d, file="res.D.add_b3.rdata")
```

Final plots
``` {r}
# 1. set env ----
library(tidyverse)
library(ggtext)
setwd("/Users/sunxin/Documents/Projects/2019_Tiger_ancient_historical/South_China_tiger/Progress_project/3.Gene_flow/2023_adm2/")

# Plot functioon ----
plot_mig_group = 
    function(list_plot, pop_fix, res_d,
             pop_order, shape_pop,
             color_pop="ppX",
             plot_d=FALSE, pfix=TRUE,
             do_plot=TRUE){
        # input:
        #   list_plot:  list of pop str to plot
        #   pop_fix:    df of pop with fixed ind
        #   res_d:      d stat result
        #   pop_order:  plot pop_order
        #   shape_pop:  pop to plot with shape, pp2, pp3
        #   color_pop:  pop to plot with color, default=popX
        #   plot_d:     plot d or z, default z
        #   pfix:       fix p2-3
        #   do_plot:    if plot, default=TRUE
        # output:
        #   ggplot obj
        
        # browser()
        df_plot = c()
        for(i in 1:length(list_plot)){
            tmp_run_g = list_plot[i]
            if (res_d %>% filter(run_g==tmp_run_g) %>% count() %>% pull(n) == 0){
                # use flip pop3, pop4
                tmp_res = res_d %>%
                    filter(run_g_f==tmp_run_g) %>%
                    mutate(est=-est,
                           z=-z,
                           popX=pop3,
                           pop3=pop4,
                           pop4=popX,
                           ppX=pp3,
                           pp3=pp4,
                           pp4=ppX,
                           run_g_p=tmp_run_g)
            }else{
                tmp_res = res_d %>%
                    filter(run_g==tmp_run_g) %>%
                    mutate(popX=pop4,
                           ppX=pp4,
                           run_g_p=tmp_run_g)
            }
            if(isTRUE(pfix)){
                # fix p2-3
                tmp_res = tmp_res %>%
                    filter(pop2==pop_fix[unlist(str_split(tmp_run_g, ","))[2],],
                           pop3==pop_fix[unlist(str_split(tmp_run_g, ","))[3],],
                    )
            }
            tmp_res = tmp_res %>%
                mutate(is_sig = case_when(abs(z)>=3 ~ "Yes",
                                          TRUE ~ "No"))
            df_plot = rbind(df_plot, tmp_res)
        }
        
        p = df_plot 
        if (isFALSE(plot_d)){
            # plot z
        }else{
            # plot d
            if(isTRUE(do_plot)){
                p = df_plot %>%
                    filter(pp1_2=="PUN") %>%
                    group_by(est) %>%
                    arrange(est) %>%
                    group_by(ppX) %>%
                    mutate(ppX = factor(ppX, levels=pop_order)) %>%
                    mutate(popX = factor(popX, levels=unique(popX))) %>% 
                    
                    ggplot() +
                    geom_pointrange(
                        aes(x=popX,y=est,
                            ymin=est-se, ymax=est+se,
                            fill=.data[[color_pop]],
                            shape=.data[[shape_pop]],
                            alpha=is_sig),
                        size=0.8, stroke=0.2,
                        position=position_dodge(width = 1),
                        linetype=1) +
                    geom_hline(yintercept = 0, linetype = 2) 
                
            }
        }
        return(p)
    }



# 2. final plot ----
## 2.0 load data ----

# load meta, sample meta, legend meta
# load plot function

# load data
load("./res_d.rdata")
# load("./res.D.add_b1.rdata")
# load("./res.D.add_b2.rdata")
# load("./res.D.add_b3.rdata")

res_d[res_d[,]=="pti217"] = "Pti217"
res_d[res_d[,]=="PtiP121"] = "Ptip121"

res_d_g = res_d %>%
    mutate(pp1 = "OUT",
           pp1_2 = case_when(pop1 == "SRR836361" ~ "PLE",
                             pop1 == "SRR836372" ~ "PUN"),
           pp2 = sample_meta[pop2,]$pop,
           pp3 = sample_meta[pop3,]$pop,
           pp4 = sample_meta[pop4,]$pop,
    ) %>%
    mutate(run_g = str_c(pp1,pp2,pp3,pp4, sep = ","),
           run_g_f = str_c(pp1,pp2,pp4,pp3, sep = ","),
    ) %>% 
    filter(! is.na(est))

pop_fix = sample_meta %>%
    filter(d_fix=="Yes") %>%
    select(c(id,pop)) %>%
    column_to_rownames("pop")


## 2.1 Fig. 4A ----
# D(Out, TIG, VIR, PopX)
list_mig3 = c('OUT,TIG,VIR,COR',
              'OUT,TIG,VIR,JAX',
              'OUT,TIG,VIR,AMO',
              'OUT,TIG,VIR,ALT',
              'OUT,TIG,VIR,RUSA')

p = plot_mig_group(
    list_plot = list_mig3,
    pop_fix=pop_fix,
    res_d=res_d_g,
    pop_order = pop_order,
    shape_pop = "pp2",
    pfix = FALSE,
    plot_d=FALSE)

color_pop="ppX"
shape_pop = "pp4"
p %>% 
    # fix pp2,3
    filter(pop2=="pti331") %>% 
    filter(pop3=="PTV06") %>%
    mutate(popX=pop4,
           ppX=pp4) %>% 
    
    filter(pp1_2=="PUN") %>%
    group_by(est) %>%
    arrange(est) %>%
    group_by(ppX) %>%
    mutate(ppX = factor(ppX, levels=pop_order)) %>%
    mutate(popX = factor(popX, levels=unique(popX))) %>% 
    
    # view()
    ggplot() +
    geom_pointrange(
        aes(x=ppX,y=est,
            ymin=est-se, ymax=est+se,
            fill=is_sig),
        size=0.8,
        shape=21,
        stroke=0.2,
        position=position_jitter(width = 0.2),
        linetype=1) +
    geom_hline(yintercept = 0, linetype = 2) +
    
    p_theme_d +
    
    # scale_fill_manual(values = c("#999999", "white"),
    # breaks= c("Yes", "No")) +
    scale_fill_manual(values = c("white", "#999999"),
                      breaks= c("No","Yes")) +
    
    theme(legend.position = "none",
          # legend.box = "vertical",
          aspect.ratio = 3/4) + 
    coord_flip() +
    guides(
        fill=guide_legend(
            title="PopX", override.aes = list(shape=21)),
        shape=guide_legend(
            title="Pop3", override.aes = list(fill="black")),
        alpha=guide_legend(
            title="zscore > 3"))

## 2.2 Fig. 4B ----
# D(OUT, PopX, ALT, RUSA)
list_mig3 = c('OUT,SUM,ALT,RUSA',
              'OUT,TIG,ALT,RUSA',
              'OUT,COR,ALT,RUSA',
              'OUT,JAX,ALT,RUSA',
              'OUT,VIR,ALT,RUSA'
              )

p = plot_mig_group(
    list_plot = list_mig3,
    pop_fix=pop_fix,
    res_d=res_d_g,
    pop_order = pop_order,
    shape_pop = "pp2",
    pfix = FALSE,
    plot_d=FALSE)

color_pop="ppX"
shape_pop = "pp2"
p %>% 
    # fix pp2,3
    filter(pop3 %in% c("RFET0003")) %>%
    filter(pop4=="RUSA21") %>%
    mutate(popX=pop2,
           ppX=pp2) %>% 
    
    filter(pp1_2=="PUN") %>%
    group_by(est) %>%
    arrange(est) %>%
    group_by(ppX) %>%
    mutate(ppX = factor(ppX, levels=pop_order)) %>%
    mutate(popX = factor(popX, levels=unique(popX))) %>% 
    
    # view()
    ggplot() +
    geom_pointrange(
        aes(x=ppX,y=est,
            ymin=est-se, ymax=est+se,
            group=pp3,
            shape=pp3,
            fill=is_sig),
        size=0.8,
        stroke=0.2,
        position=position_jitter(width = 0.4, seed = 123),
        # position=position_jitterdodge(
        #     dodge.width = 1),
        linetype=1) +
    geom_hline(yintercept = 0, linetype = 2) +
    
    p_theme_d +
    
    # scale_fill_manual(values = c("#999999", "white"),
    # breaks= c("Yes", "No")) +
    scale_fill_manual(values = c("white", "#999999"),
                      breaks= c("No","Yes")) +
    
    scale_shape_manual(values = meta_legend$shape,
                       breaks = meta_legend$pop) +
    
    theme(legend.position = "none",
          legend.title = element_markdown(),
          aspect.ratio = 3/4) + 
    coord_flip() +
    guides(
        fill=guide_legend(
            title="*D*-statistics is significant"),
    )

ggsave("D.mig_rusa_rev2_alt-rusa.grey.6x9.pdf", width = 6, height = 9, device = "pdf", dpi = 500)
ggsave("D.mig_rusa_rev2_alt-rusa.grey.legend.6x9.pdf", width = 6, height = 9, device = "pdf", dpi = 500)

## 2.3 Fig. S5A ----
# D(X, AMO_1, AMO_X)
# load add2
list_mig4 = c(
    'OUT,COR,AMO,AMO',
    'OUT,JAX,AMO,AMO',
    'OUT,SUM,AMO,AMO',
    'OUT,TIG,AMO,AMO')

p = plot_mig_group(
    list_plot = list_mig4,
    pop_fix=FALSE,
    pfix=FALSE,
    res_d=res_d_g,
    pop_order = pop_order,
    shape_pop = "pp2",
    color_pop = "pp2",
    do_plot=FALSE)


popX_order = p %>% 
    # order by the mean est of popX
    group_by(popX) %>%
    summarise(m=mean(est)) %>% 
    arrange(m) %>% 
    pull(popX)

p %>%
    filter(pp1_2=="PUN") %>%
    group_by(est) %>%
    arrange(est) %>%
    group_by(ppX) %>%
    mutate(ppX = factor(ppX, levels=pop_order)) %>%
    mutate(popX = factor(popX, levels=popX_order)) %>% 
    
    ggplot() +
    geom_pointrange(
        aes(x=popX,y=est,
            ymin=est-se, ymax=est+se,
            fill=pp2,
            alpha=is_sig),
        shape=21,
        size=0.8, stroke=0.2,
        position=position_dodge(width = 1),
        linetype=1) +
    geom_hline(yintercept = 0, linetype = 2) +
    
    p_theme_d +
    
    scale_fill_manual(values = meta_legend$color,
                      breaks = meta_legend$pop) +
    scale_shape_manual(values = meta_legend$shape,
                       breaks = meta_legend$pop) +
    scale_alpha_discrete(range = c(0.3, 1),
                         breaks= c("No", "Yes")) +
    
    theme(legend.position = "right",
          aspect.ratio = 2) + 
    coord_flip() +
    guides(
        fill=guide_legend(
            title="Pop2", override.aes = list(shape=21)),
        # shape=guide_legend(
            # title="Pop3", override.aes = list(fill="black")),
        alpha=guide_legend(
            title="zscore > 3", override.aes = list(fill="black")))


ggsave("D.X_amo_amo.6x9.png", width = 6, height = 9, device = "png", dpi = 500)
ggsave("D.X_amo_amo.6x9.pdf", width = 6, height = 9, device = "pdf", dpi = 500)

## 2.4 Fig. S5B ----
# load data
list_mig1 = c('OUT,SUM,COR,TIG',
              'OUT,SUM,COR,JAX',
              'OUT,SUM,JAX,COR',
              'OUT,SUM,JAX,TIG',
              'OUT,SUM,COR,AMO',
              'OUT,SUM,JAX,AMO',
              'OUT,SUM,COR,ALT',
              'OUT,SUM,JAX,ALT',
              'OUT,SUM,COR,VIR',
              'OUT,SUM,JAX,VIR',
              'OUT,SUM,COR,RUSA',
              'OUT,SUM,JAX,RUSA'
)


p = plot_mig_group(
    list_plot = list_mig1,
    pop_fix=pop_fix,
    res_d=res_d_g,
    pop_order = pop_order,
    shape_pop = "pp3",
    plot_d=TRUE)

p + p_theme_d +
    scale_fill_manual(values = meta_legend$color,
                      breaks = meta_legend$pop) +
    scale_shape_manual(values = meta_legend$shape,
                       breaks = meta_legend$pop) +
    scale_alpha_discrete(range = c(0.3, 1),
                         breaks= c("No", "Yes")) +
    
    theme(legend.position = "right",
          aspect.ratio = 2) + 
    coord_flip() +
    guides(
        fill=guide_legend(
            title="PopX", override.aes = list(shape=21)),
        shape=guide_legend(
            title="Pop3", override.aes = list(fill="black")),
        alpha=guide_legend(
            title="zscore > 3"))

ggsave("D.mig_sum.6x9.png", width = 6, height = 9, device = "png", dpi = 500)





## 2.5 Fig. S5C ----
list_mig1 = c('OUT,SUM,TIG,COR',
              'OUT,SUM,TIG,JAX',
              'OUT,SUM,TIG,AMO',
              'OUT,SUM,TIG,AMUR',
              'OUT,SUM,TIG,VIR',
              'OUT,SUM,TIG,RUSA'
)


p = plot_mig_group(
    list_plot = list_mig1,
    pop_fix=pop_fix,
    res_d=res_d_g,
    pop_order = pop_order,
    shape_pop = "pp3",
    plot_d=FALSE)

p %>%
    
    filter(pp1_2=="PUN") %>%
    group_by(est) %>%
    arrange(est) %>%
    group_by(ppX) %>%
    mutate(ppX = factor(ppX, levels=pop_order)) %>%
    mutate(popX = factor(popX, levels=unique(popX))) %>% 
    
    ggplot() +
    geom_pointrange(
        aes(x=popX,y=est,
            ymin=est-se, ymax=est+se,
            fill=ppX,
            alpha=is_sig),
        size=0.8, stroke=0.2,
        shape=21,
        position=position_dodge(width = 1),
        linetype=1) +
    
    geom_hline(yintercept = 0, linetype = 2) +
    
    p_theme_d +
    scale_fill_manual(values = meta_legend$color,
                      breaks = meta_legend$pop) +
    scale_alpha_discrete(range = c(0.3, 1),
                         breaks= c("No", "Yes")) +
    
    theme(legend.position = "right",
          aspect.ratio = 2) + 
    coord_flip() +
    guides(
        fill=guide_legend(
            title="PopX", override.aes = list(shape=21)),
        alpha=guide_legend(
            title="zscore > 3", override.aes = list(fill="black")))

ggsave("D.mig_sum_tig.6x9.png", width = 6, height = 9, device = "png", dpi = 500)

```

## 2. AdmixtureGraph
AdmixtureGraphs were explored with qpBrute.

par file `par_qpg`
```
DIR:            .
SSS:            auto_neutral_withOUT_hq_trvn
genotypename:   DIR/SSS.geno
snpname:        DIR/SSS.snp_dis
indivname:      DIR/SSS.indi_indi
outpop:         NULL
blgsize:        0.005
lsqmode:        YES
diag:           .0001
hires:          YES
```

The base graph used is (PUN,(PLE, SUM)) `out2_sum-0d813314d37b.graph`.  
```
qpBrute \
    --par  par_qpg\
    --prefix out2_sum2 \
    --pops PLE,PUN,SUM1\
    --out out2_sum

```

The remaining pops [COR1', 'JAX1', 'RUSA', 'TIG1', 'AMO1', 'PTV1', 'AMUR1'] were added to the graph with qpBrute with default settings.
```
qpBrute \
    --par  par_qpg\
    --prefix out2_sum2 \
    --pops COR1,JAX1,RUSA,TIG1,AMO1,PTV1,AMUR1 \
    --out Out \
    --qpgraph out2_sum-0d813314d37b.graph
```
Only one fitted graph left after the search.