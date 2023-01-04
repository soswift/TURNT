# process each run using R script
setwd("~/Documents/Bioinformatics/Projects/TURNT/data/fdom_data/")
source("../../src/process_aqualog_functions.R")



# Process individual runs -------------
# round 1
run1 <- process_aqualog(data_directory ="matrices",
                        run_name = "TURNT",
                        sample_key_file = "TURNT_fdom_run_20220714.tsv")

saveRDS(run1, "run1.rds")

#  Plot TURNT eems -------------------

  

