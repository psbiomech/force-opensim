### FORCE - SDP OPENSIM DESCRIPTIVES BY SUBJECT
#
# Prasanna Sritharan, May 2022


library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/data/FORCe/outputdatabase/csvfolder"
outfolder <- "r-output"

# load OpenSim results
osimdatafile <- "force_sdp_results_all_trials.csv"
osim  <- read_csv(file.path(srcfolder, osimdatafile))


# ******************************
# CLEAN RAW DATA


# bind subject info to OpenSim results, recode/rename/rearrange for readability
osim <- osim %>% 
          mutate(group=if_else(grepl("C", subject), "con", "frc"), .after=subject) %>% 
          mutate(aff_side=recode(aff_side, `0`="C", `1`="R", `2`="L", `3`="B", `-1`="N"), 
                 dom_foot=recode(dom_foot, `1`="R", `2`="L"),
                 data_leg=toupper(data_leg)) %>% 
          mutate(sex=recode(sex, `1`="M", `2`="F"), .after=sex)

# determine shomri score for the data leg
osim <- osim %>%
          mutate(shomri_data_leg=if_else(data_leg=="R", shomri_r, shomri_l), .after=shomri_l)

# version 1: determine type of leg for the data: dom, ndom, sym or asym
osim <- osim %>%
  mutate(data_leg_type=if_else((data_leg==dom_foot) & (aff_side=="C"), "dom",
                       if_else((data_leg!=dom_foot) & (aff_side=="C"), "ndom",
                       if_else((data_leg==aff_side) | (aff_side=="B"), "sym", "asym"))),
                       .after=data_leg)

# version 2: determine type of leg for the data: ctrl, more or less
# Assume: if no shomri data, dominant foot is the more affected limb 
osim <- osim %>%
  mutate(data_leg_type2=if_else((aff_side=="C"), "ctrl",
                        if_else((data_leg==aff_side), "more", 
                        if_else((aff_side=="B") & (data_leg=="R") & (shomri_r>=shomri_l), "more", 
                        if_else((aff_side=="B") & (data_leg=="L") & (shomri_l>=shomri_r), "more", 
                        if_else((aff_side=="B") & ((shomri_r=="NA" | shomri_l=="NA")) & (dom_foot=="R") & (data_leg=="R"), "more", 
                        if_else((aff_side=="B") & ((shomri_r=="NA" | shomri_l=="NA")) & (dom_foot=="L") & (data_leg=="L"), "more", "less")))))),
                        .after=data_leg_type)


# ******************************
# NORMALISE OPENSIM RESULTS

# normalisation factors
g = 9.81
osimnorm <- osim %>% 
              mutate(bw=mass*g, bwht=bw*height, .after=height)


# joint moments
osimnorm <- osimnorm %>% 
              mutate(across(t1:t101, ~ if_else(analysis=="id" & grepl("pelvis_t.", variable), as.numeric(.x)/bw, 
                                       if_else(analysis=="id" & variable!="time", as.numeric(.x)*100/bwht, as.numeric(.x)))))

# muscle forces
osimnorm <- osimnorm %>% 
              mutate(across(t1:t101, ~ if_else(analysis=="so" & variable!="time", as.numeric(.x)/bw, as.numeric(.x))))

# joint forces
osimnorm <- osimnorm %>% 
              mutate(across(t1:t101, ~ if_else(analysis=="jr" & grepl(".+_f[xyz]$", variable), as.numeric(.x)/bw, 
                                       if_else(analysis=="jr" & grepl(".+_m[xyz]$", variable), as.numeric(.x)*100/bwht, as.numeric(.x)))))
                     

                     
# ******************************
# CALCULATE DESCRIPTIVES PER SUBJECT

# group by subject and all non-numeric fields, count trials in each group
osimnorm <- osimnorm %>% 
              group_by(subject, group, sex, dom_foot, aff_side, task, data_leg, data_leg_type, data_leg_type2, data_leg_role, analysis, variable) %>% 
              mutate(ntrials=n(), .before=analysis)

# descriptives for temporal variables
subj_mean <- osimnorm %>% 
                summarise_if(is.numeric, mean) %>% 
                mutate(statistic="mean", row_number=row_number())
subj_sd <- osimnorm %>% 
                summarise_if(is.numeric, sd) %>% 
                mutate(statistic="sd", row_number=row_number())

# bind and interleave mean and sd rows
descriptives <- subj_mean %>% 
                  bind_rows(subj_sd) %>% 
                  arrange(row_number, .by_group=TRUE) %>% 
                  select(-row_number)

# rearrange columns for readability
descriptives <- descriptives %>% 
                  relocate(c(age, mass, height), .before=sex) %>% 
                  relocate(task, .after=subject) %>% 
                  relocate(c(ntrials, analysis, variable, statistic), .before=t1) %>% 
                  relocate(c(shomri_r, shomri_l), .after=sex) %>% 
                  relocate(shomri_data_leg, .after=data_leg_role) %>% 
                  select(-c(bw, bwht))


# write to file
outpath = file.path(srcfolder, outfolder)
if (!dir.exists(outpath)) {dir.create(outpath)}
write_csv(descriptives, file.path(outpath, "force_sdp_results_descriptives_by_subject.csv"))


