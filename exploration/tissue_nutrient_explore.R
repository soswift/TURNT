# exploration of nutrient data for turbinaria tissue and water samples
#Tissue N and C -----
# data received from Popp lab 05/23/22
 library(data.table)
 library(ggplot2)
 library(ggpubr)
 library(rstatix)
 library(GGally)

 setwd("~/Documents/Bioinformatics/Projects/TURNT/exploration/")
 
  # read
  t_meta_file = "../data/nutrient_data/tissue/TURNT_sample_sheet - tissue_nutrient_20220526.csv"
  raw_nut_file = "../data/nutrient_data/tissue/turnt_tissue_nutrient_20220523.csv"
 
  sgd_nut_file = "../data/nutrient_data/tissue/sgd_tissue.csv"
  
  t_meta = fread(t_meta_file)
  t_dat = fread(raw_nut_file)
 
  sgd_dat = fread(sgd_nut_file)
  
  
  
  
  # organize
  t_nut = merge(t_meta, t_dat, by = "individual")
  t_nut[, perc_N :=  ug_N/(10* weight_mg)]
  t_nut[, perc_C :=  ug_C/(10* weight_mg)]
  t_nut[, C_to_N := ug_C/ug_N]
  
  t_nut[, clone := gsub("(..)-[12]", "\\1",individual)]
  t_nut[, size := ifelse(grepl("-1", individual),
                         "large",
                         "small")]
  

  p_vars = c(
    "weight_mg",
    "perc_C",
    "C13_perc_vs_vpdb",
    "perc_N",
    "N_15_perc_vs_air",
    "C_to_N"
  )
  
  turnt_clrs = c(
    control = "#73d9f0",
    nutrient = "#f7b301",
    recharge = "#eb773e",
    blank = "gray"
  )


  # plot stat ----------------------------
  ggpairs(t_nut, columns = 6:12, aes(color = treatment, fill = treatment))+
    scale_color_manual(values = turnt_clrs)+
    scale_fill_manual(values = turnt_clrs)
    
  ggsave("plots/nutrient/turnt_var_plot.pdf")
  
  ggpairs(t_nut, columns = 6:12)
  
  ggsave("plots/nutrient/turnt_all_var_plot.pdf")
  
  ggpairs(sgd_dat, columns = 3:8)
  
  ggsave("plots/nutrient/sgd_var_plot.pdf")
  
  t_p_list = lapply(p_vars, function(pvar) {
      # Tukey
      stat.test = aov(as.formula(paste0(pvar, "~", "treatment")), data = t_nut)
      stat.test = tukey_hsd(stat.test)
      
      # plot
      y_max = max(t_nut[[pvar]])
      print(y_max)
      
     p= ggboxplot(t_nut,
                  x = "treatment",
                  y = pvar,
                  fill = "treatment",
                  outlier.shape = NA)+
       geom_jitter(width = 0.2, pch = 21)+

      stat_pvalue_manual(stat.test, label = "p.adj", 
                         y.position = c(y_max +1,
                                        y_max + 1.75,
                                        y_max + 2.5),
                         size = 3)+
       scale_fill_manual(values = turnt_clrs)+
       theme(legend.position = "none")+
        xlab("")
     ggpar(p, ylim = c(min(t_nut[[pvar]]), y_max + 3))
  }
  )
  
  names(t_p_list) = p_vars
  
  ggarrange(plotlist = t_p_list, ncol = 1)
  
  
  # all plots arranged
  ggsave(filename = "plots/nutrient/tissue_nutrient_by_treatment_stat.png",
         height = 10, width = 6)
  
  
 # just C:N ratio
  ggsave(filename = "plots/nutrient/C_to_N_box.png",
         plot = t_p_list$C_to_N,
         height = 7, width = 7)
  

# model stat
  t_anova = lapply(p_vars, function(p_var) {
    mod = lm(as.formula(paste0(p_var, "~ treatment")), data = t_nut)
    anova(mod)
  })
  
  names(t_anova) = p_vars
  t_anova$C_to_N
  
# Alternative units of analysis
# mean of clone  --------------------
  m_t_nut = t_nut[ , lapply(.SD, function(x) ifelse(is.character(x),
                                                    unique(x),
                                                    mean(x))),
                   .SDcols = c("treatment", p_vars),
                   by = clone]
  
 # pair plot
  ggpairs(m_t_nut, columns = 3:8, aes(color = treatment))+      
    scale_color_manual(values = turnt_clrs)

  ggsave(filename = "plots/nutrient/mean_tissue_nutrient_pair_plot.png",
         width = 10,
         height = 10)
  
# box plot
  
  m_t_p_list = lapply(p_vars, function(pvar) {
    # Tukey
    stat.test = aov(as.formula(paste0(pvar, "~", "treatment")), data = m_t_nut)
    stat.test = tukey_hsd(stat.test)
    
    # plot
    var_dat = m_t_nut[[pvar]]
    
    y_max = max(var_dat)
    y_mean = mean(var_dat)
    y_min = min(var_dat)

    
    p= ggboxplot(m_t_nut,
                 x = "treatment",
                 y = pvar,
                 fill = "treatment",
                 outlier.shape = NA)+
      geom_jitter(pch = 21, width = 0.1)+
      
      stat_pvalue_manual(stat.test, label = "p.adj", 
                         y.position = c(y_max + abs(0.1*y_mean),
                                        y_max + abs(0.2*y_mean),
                                        y_max + abs(0.3*y_mean)),
                         size = 3)+
      
      scale_fill_manual(values = turnt_clrs)+
      
      theme(legend.position = "none")+
      xlab("")
    ggpar(p, ylim = c(y_min, y_max + abs(0.4*y_mean)))
  }
  )
  
  names(m_t_p_list) = p_vars
  
  ggarrange(plotlist = m_t_p_list, ncol = 2, nrow =3)
  ggsave(filename = "plots/nutrient/mean_tissue_nutrient_by_treatment.png",
         bg = "white",
         width = 10, height = 10)

  for(i in p_vars){
   ggsave(paste0("plots/nutrient/box_",i,".png"),
          plot = m_t_p_list[[i]],
          width = 7, height = 5,
          scale = 0.75)
  }
    

# write out data frame
  fwrite(t_nut, "../data/nutrient_data/clean_tissue_nutrients_TURNT.csv")
  
  