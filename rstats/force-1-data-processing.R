### FORCE - OPENSIM PROCESSING
#
# Prasanna Sritharan, May 2022


library(tidyverse)


# folders
srcfolder <- "C:/Users/Owner/Documents/data/FORCe/outputdatabase/csvfolder"
outfolder <- "r-output"

osimdatafile <- "force_sdp_results_all_sdp.csv"
osim  <- read_csv(file.path(srcfolder, osimdatafile))

# add trial leg, and recode sex to numeric
osim <- osim %>% 
          mutate(first_event=str_extract(osim[["events_labels"]], "\\wFO"), step_leg=substr(first_event, 1, 1), .before=foot) %>% 
          select(-c("first_event", "condition"))

# load participant data spreadsheet
subjdatafile <- "FORCE-ParticipantData-All.txt"
subjdata <- read_delim(file.path(srcfolder, subjdatafile), delim="\t")

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

# bind to osim data, rearrange and recode affected and dom_foot, for sex create
# new variable as the numeric codes are required for summarise
osim <- osim %>% 
          select(-c("leg_task")) %>% 
          bind_cols(subjinfomat) %>% 
          relocate(names(subjdata)[-1], .after=group) %>% 
          mutate(group=if_else(grepl("CRT",subject), "CON", "FRC"),
                 aff_side=recode(aff_side, `0`="C", `1`="R", `2`="L", `3`="B", `-1`="N"), 
                 dom_foot=recode(dom_foot, `1`="R", `2`="L"),
                 foot=toupper(foot)) %>% 
          mutate(sex0=recode(sex, `1`="M", `2`="F"), .after=sex)

# label trials better
osim <- osim %>% 
          mutate(data_limb=if_else(tolower(step_leg)==tolower(foot), "IPSI", "CONTRA"),
                 ipsi_limb=if_else(tolower(aff_side)==tolower(step_leg) | (aff_side=="B"), "SYM", 
                            if_else((group=="CON") & tolower(dom_foot)==tolower(step_leg), "DOM", 
                            if_else((group=="CON") & tolower(dom_foot)!=tolower(step_leg), "NDOM", "ASYM"))),
                 contra_limb=if_else(ipsi_limb=="ASYM", "SYM", 
                             if_else(ipsi_limb=="DOM", "NDOM",
                             if_else(ipsi_limb=="NDOM", "DOM", "ASYM")))) %>% 
          relocate(c("ipsi_limb", "contra_limb", "data_limb"), .after=task)

# convert event times to steps, and bind
events <- osim %>% select(events_times)
eventmat <- NULL
for (ev in 1:nrow(events)) {
  erow <- as.numeric(as_vector(str_split(str_sub(events[[1]][ev],2,-2), "; ")))
  estep <- round(101 * (erow - erow[1]) / (erow[6] - erow[1]))
  eventmat <- rbind(eventmat, estep)
}
eventmat <- as_tibble(eventmat)
colnames(eventmat) <- lapply(1:6, function(x) {paste0("E", x)})
osim <- osim %>% 
          bind_cols(eventmat) %>% 
          relocate(colnames(eventmat), .before=analysis)

# write to file
write_csv(osim, file.path(srcfolder, "force_sdp_results_updated_sdp.csv"))


