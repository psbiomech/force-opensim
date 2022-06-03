### FORCE - DESCRIPTIVES
#
# Prasanna Sritharan, May 2022

library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/projects/force-moco/rstats"
datafolder <- "r-data"
outfolder <- "r-output"

# updated OpenSim data spreadsheet
osimdatafile <- "force_sdp_results_normalised_sdp.csv"
osim  <- read_csv(file.path(srcfolder, datafolder, osimdatafile))


# ******************************
# WITHIN SUBJECT MEAN AND SD

# group
osim <- osim %>% group_by(subject, data_limb, ipsi_limb, contra_limb, analysis, variable)

# descriptives for temporal variables
summary_mean_ws <- osim %>% 
                  summarise_if(is_numeric, mean) %>% 
                  mutate(statistic="mean", row_number=row_number())
summary_sd_ws <- osim %>% 
                  summarise_if(is_numeric, sd) %>% 
                  mutate(statistic="sd", row_number=row_number())

# bind and interleave
within_subject <- summary_mean_ws %>% 
                    bind_rows(summary_sd_ws) %>% 
                    arrange(row_number, .by_group=TRUE) %>% 
                    relocate(c("age", "mass", "height"), .before=analysis) %>% 
                    relocate(statistic, .after=variable) %>% 
                    select(-row_number)



# ******************************
# BETWEEN GROUP MEAN AND SD

# group, remove sd rows
within_subject_data <- within_subject %>% 
                          ungroup() %>% 
                          filter(statistic=="mean") %>% 
                          group_by(data_limb, ipsi_limb, contra_limb, analysis, variable)

# descriptives for temporal variables
summary_mean_bg <- within_subject_data %>% 
                    summarise_if(is_numeric, mean) %>% 
                    mutate(statistic="mean", row_number=row_number())
summary_sd_bg <- within_subject_data %>% 
                    summarise_if(is_numeric, sd) %>% 
                    mutate(statistic="sd", row_number=row_number())

# bind and interleave, remove non-temporal variables
between_group <- summary_mean_bg %>% 
                    bind_rows(summary_sd_bg) %>% 
                    arrange(row_number, .by_group=TRUE) %>% 
                    relocate(c("age", "mass", "height"), .before=analysis) %>% 
                    relocate(statistic, .after=variable) %>% 
                    select(-row_number)
                    



# ******************************
# WRITE TO CSV

# within subject
within_subject <- within_subject %>% select(-c("age", "mass", "bw", "bwht"))
write_csv(within_subject, file.path(srcfolder, outfolder, "force_sdp_normalised_descriptives_subject.csv"))

# between group
between_group <- between_group %>% select(-c("age", "mass", "bw", "bwht"))
write_csv(between_group, file.path(srcfolder, outfolder, "force_sdp_normalised_descriptives_group.csv"))