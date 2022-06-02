### FORCE - DESCRIPTIVES
#
# Prasanna Sritharan, May 2022

library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/projects/force-moco/rstats"
datafolder <- "r-data"
outfolder <- "r-output"

# updated OpenSim data spreadsheet
osimdatafile <- "force_sdp_results_updated_sdp.csv"
osim  <- read_csv(file.path(srcfolder, datafolder, osimdatafile))


# group
osim <- osim %>% group_by(data_limb, ipsi_limb, contra_limb, analysis, variable)

# descriptives for temporal variables
summary_mean <- osim %>% 
                  summarise_if(is_numeric, mean) %>% 
                  mutate(statistic="mean", row_number=row_number())
summary_sd <- osim %>% 
                summarise_if(is_numeric, sd) %>% 
                mutate(statistic="sd", row_number=row_number())

# bind and interleave, remove non-temporal variables
summary_descriptives <- summary_mean %>% 
                          bind_rows(summary_sd) %>% 
                          arrange(row_number, .by_group=TRUE) %>% 
                          relocate(c("age", "mass", "height"), .before=analysis) %>% 
                          relocate(statistic, .after=variable) %>% 
                          select(-c("row_number", "age", "height", "mass"))
    
# write to file
write_csv(summary_descriptives, file.path(srcfolder, outfolder, "force_sdp_descriptives.csv"))

