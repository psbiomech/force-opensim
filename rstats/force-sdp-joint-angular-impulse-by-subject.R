### FORCE - JOINT ANGULAR IMPULSE DESCRIPTIVES
# Prasanna Sritharan, July 2022



library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/data/FORCe/outputdatabase/csvfolder"
outfolder <- "r-output"

# load data
osimdatafile <- "force_sdp_joint_angular_impulse_all.csv"
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
          select(-c(`group`, `id`)) %>% 
          mutate(sex=recode(sex, `1`="M", `2`="F"), .after=sex) %>%
          mutate(group=if_else(grepl("CRT", subject), "CON", "FRC"), .after=subject) %>% 
          mutate(aff_side=recode(aff_side, `0`="C", `1`="R", `2`="L", `3`="B", `-1`="N"), 
            dom_foot=recode(dom_foot, `1`="R", `2`="L"),
            foot=toupper(foot),
            `type`=recode(type, `net`="full", `windows`="window"),
            `window`=recode(window, `net`="full")) %>% 
          rename(`period`=`type`)


# normalise data
grav = 9.81
osimnorm <- osim %>% mutate(net=(100*net)/(mass*grav*height),
                            pos=(100*pos)/(mass*grav*height),
                            neg=(100*neg)/(mass*grav*height))



# ******************************
# CALCULATE DESCRIPTIVES





          
