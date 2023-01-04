# exploration of nutrient data for turbinaria tissue and water samples
#Tissue N and C -----
# data received from Popp lab 05/23/22
 library(data.table)
 library(ggplot2)
 library(cowplot)
 

 setwd("~/Documents/Bioinformatics/Projects/TURNT/exploration/")
 
  # read
  t_meta_file = "../data/nutrient_data/TURNT_sample_sheet - tissue_nutrient_20220524.csv"
  raw_nut_file = "../data/nutrient_data/turnt_tissue_nutrient_20220523.csv"
 
  t_meta = fread(t_meta_file)
  t_dat = fread(t_nut_file)
 
  # organize
  t_nut = merge(t_meta, t_dat, by = "individual")
  t_nut[, perc_N := ug_N/weight_mg]
  t_nut[, perc_C := ug_C/weight_mg]
 
  # plot
  p_vars = c(
    "weight_mg",
    "C13_perc_vs_vpdb",
    "N_15_perc_vs_air",
    "perc_N"
  )
 
  t_p_list = lapply(p_vars, function(p_var) {
   ggplot(t_nut, aes_string(x = "treatment",
                            y = p_var,
                            fill = "treatment"))+
    geom_boxplot()+
    geom_jitter(width = 0.1,
                shape = 21)+
      theme_minimal()+
      xlab("")+
      theme(legend.position = 'none')
  }
 )

  
  t_p = plot_grid(plotlist = t_p_list,ncol = 1)
  t_p
  ggsave(plot = t_p, filename = "plots/tissue_nutrient_by_treatment.png", bg = "white")
  

