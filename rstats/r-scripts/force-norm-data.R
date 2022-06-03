### FORCE - NORMALISE DATA
#
# Prasanna Sritharan, June 2022

library(tidyverse)

# folders
srcfolder <- "C:/Users/Owner/Documents/projects/force-moco/rstats"
datafolder <- "r-data"
outfolder <- "r-output"


# load OpenSim data spreadsheet
osimdatafile <- "force_sdp_results_updated_sdp.csv"
osim  <- read_csv(file.path(srcfolder, datafolder, osimdatafile))


# normalisation factors
g = 9.81
osimnorm <- osim %>% 
              mutate(bw=mass*g, bwht=bw*height, .after=height)

# parse data line by line and normalise if require (this is inefficient, need
# to find a way to use tidyverse to do this quickly)
tcolidx = which(grepl("^t\\d+", names(osimnorm)))
rows = nrow(osimnorm)
for (r in 1:rows) {
  
  # inverse dynamics
  if ((osimnorm$analysis[r]=="id") & (osimnorm$variable[r]!="time")) {
    if (grepl("^pelvis_t\\w_force", osimnorm$variable[r])) {
      osimnorm[r, tcolidx] = osimnorm[r, tcolidx] / osimnorm$bw[r] 
    } else {
      osimnorm[r, tcolidx] = 100 * osimnorm[r, tcolidx] / osimnorm$bwht[r]  
    }
    
  }

}
    

# write to file
write_csv(osimnorm, file.path(srcfolder, datafolder, "force_sdp_results_normalised_sdp.csv"))
