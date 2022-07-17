### FORCE - JOINT ANGULAR IMPULSE DESCRIPTIVES
# Prasanna Sritharan, July 2022



library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/data/FORCe/outputdatabase/csvfolder"
outfolder <- "r-output"

# load data
osimdatafile <- "force_sdp_joint_angular_impulse_all_trials.csv"
osim <- read_csv(file.path(srcfolder, osimdatafile))

# load participant info spreadsheet
subjdatafile <- "FORCE-ParticipantData-All.txt"
subjdata <- read_delim(file.path(srcfolder, subjdatafile), delim="\t")



# ******************************
# PREPARE DATA

# iterate and add participant info
subjinfomat <- NULL
for (r in 1:nrow(subjdata)){
  
  # get the data rows for the subject
  subj <- subjdata[["id"]][r]
  subjrows <- osim %>% filter(subject == subj) 
  
  # get the subject info and repeat rows
  subjinfo <- subjdata[r,]
  subjinforep <- subjinfo[rep(1, each=nrow(subjrows)),]
  
  # bind
  subjinfomat <- bind_rows(subjinfomat, subjinforep)
  
}

# bind to osim data, clean up and recode (mainly for readability)
osim <- osim %>% 
          bind_cols(subjinfomat) %>%
          relocate(names(subjinfomat), .after=group) %>%
          relocate(movement, .after=aff_side) %>% 
          select(-c(`group`, `id`)) %>% 
          mutate(sex=recode(sex, `1`="M", `2`="F"), .after=sex) %>%
          mutate(group=if_else(grepl("CRT", subject), "CON", "FRC"), .after=subject) %>% 
          mutate(aff_side=recode(aff_side, `0`="C", `1`="R", `2`="L", `3`="B", `-1`="N"), 
            dom_foot=recode(dom_foot, `1`="R", `2`="L"),
            data_leg=toupper(data_leg),
            period=recode(period, `net`="full", `windows`="window"),
            `window`=recode(window, `net`="full")) %>% 
          mutate(variable=str_sub(variable, 1, -8)) %>% 
          rename(dom_leg=dom_foot, task=movement)

# determine type of leg for the data: dom, ndom, sym or asym
osim <- osim %>%
  mutate(data_leg_type=if_else((data_leg==dom_leg) & (aff_side=="C"), "dom",
                       if_else((data_leg!=dom_leg) & (aff_side=="C"), "ndom",
                       if_else((data_leg==aff_side) | (aff_side=="B"), "sym", "asym"))),
                       .after=data_leg)

# determine shomri score for the data leg
osim <- osim %>%
  mutate(data_leg_shomri=if_else(data_leg=="R", r_shomri_total, l_shomri_total), .after=data_leg_role)



# ******************************
# NORMALISE OPENSIM RESULTS

# normalise data
grav = 9.81
osimnorm <- osim %>% mutate(net=(100*net)/(mass*grav*height),
                            positive=(100*positive)/(mass*grav*height),
                            negative=(100*negative)/(mass*grav*height))



# ******************************
# CALCULATE DESCRIPTIVES PER SUBJECT

# group by subject and all non-numeric fields, count trials in each group
osimnorm <- osimnorm %>% 
              group_by(subject, group, sex, dom_leg, aff_side, task, data_leg, data_leg_type, data_leg_role, period, window, variable) %>% 
              mutate(ntrials=n(), .before=period)

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
  relocate(c(ntrials, period, window, variable, statistic), .before=net) %>% 
  relocate(c(r_shomri_total, l_shomri_total), .after=sex) %>% 
  relocate(data_leg_shomri, .after=data_leg_role)


# write to file
outpath = file.path(srcfolder, outfolder)
if (!dir.exists(outpath)) {dir.create(outpath)}
write_csv(descriptives, file.path(outpath, "force_sdp_joint_angular_impulse_descriptives_by_subject.csv"))

