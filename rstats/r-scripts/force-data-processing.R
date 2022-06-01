### FORCE OPENSIM PROCESSING
#
# Prasanna Sritharan, May 2022


library(tidyverse)

# folders
srcfolder <- "C:/Users/Owner/Documents/projects/force-moco/rstats"
datafolder <- "r-data"
outfolder <- "r-output"

# load OpenSim data spreadsheet
osimdatafile <- "force_sdp_results_all_sdp.csv"
osim  <- read_csv(file.path(srcfolder, datafolder, osimdatafile))

# add trial leg
osim <- osim %>% 
          mutate(first_event=str_extract(osim[["events_labels"]], "\\wFO"), step_leg=substr(first_event, 1, 1), .before=foot) %>% 
          select(-c("first_event"))

# load participant data spreadsheet
subjdatafile <- "FORCE-ParticipantData-All.txt"
subjdata <- read_delim(file.path(srcfolder, datafolder, subjdatafile), delim="\t")



# iterate and add participant info
subjinfomat <- NULL
for (r in 1:nrow(subjdata)){
  
  # get the data rows for the subject
  subj <- subjdata[["id"]][r]
  subjrows <- osim %>% filter(subject == subj) 
  
  # get the subject info and repeat rows
  subjinfo <- subjdata[r,][-1]
  subjinforep <- subjinfo[rep(1, each=nrow(subjrows)),]
  
  # bind
  subjinfomat <- bind_rows(subjinfomat, subjinforep)
  
}

# bind to osim data, set groups, rearrange and recode affected and dom_foot
osim <- osim %>% 
          select(-c("leg_task")) %>% 
          bind_cols(subjinfomat) %>% 
          relocate(names(subjdata)[-1], .after=group) %>% 
          mutate(group=if_else(grepl("CRT",subject), "CON", "SYM"),
                 aff_side=recode(aff_side, `0`="C", `1`="R", `2`="L", `3`="B", `-1`="N"), 
                 dom_foot=recode(dom_foot, `1`="R", `2`="L"),
                 foot=toupper(foot))

# overwrite condition to be ispilateral or contralateral foot
osim <- osim %>% 
          mutate(condition=if_else(tolower(step_leg)==tolower(foot), "IPSI", "CONTRA"))



