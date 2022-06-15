### FORCE - NORMALISE DATA
#
# Prasanna Sritharan, June 2022

library(tidyverse)

# folders
srcfolder <- "C:/Users/Owner/Documents/data/FORCe/outputdatabase/csvfolder"
outfolder <- "r-output"


# load OpenSim data spreadsheet
osimdatafile <- "force_sdp_results_updated_sdp.csv"
osim  <- read_csv(file.path(srcfolder, osimdatafile))


# normalisation factors
g = 9.81
osimnorm <- osim %>% 
              mutate(bw=mass*g, bwht=bw*height, .after=height)

# parse data line by line and normalise if required (this is inefficient, need
# to find a way to use tidyverse to do this quickly, prob mutate across)
tcolidx = which(grepl("^t\\d+", names(osimnorm)))
rows = nrow(osimnorm)
for (r in 1:rows) {
  
  if (osimnorm$variable[r]!="time") {
  
    # inverse dynamics
    if (osimnorm$analysis[r]=="id") {
      if (grepl("^pelvis_t\\w_force", osimnorm$variable[r])) {
        osimnorm[r, tcolidx] = osimnorm[r, tcolidx] / osimnorm$bw[r] 
      } else {
        osimnorm[r, tcolidx] = 100 * osimnorm[r, tcolidx] / osimnorm$bwht[r]  
      }
      
    # static optimisation muscle forces
    } else if (osimnorm$analysis[r]=="so") {
      osimnorm[r, tcolidx] = osimnorm[r, tcolidx] / osimnorm$bw[r]
  
    # joint reaction forces and moments
    } else if (osimnorm$analysis[r]=="jr") {
      if (grepl(".+_f[xyz]$", osimnorm$variable[r])) {
        osimnorm[r, tcolidx] = osimnorm[r, tcolidx] / osimnorm$bw[r]
      } else {
        osimnorm[r, tcolidx] = 100 * osimnorm[r, tcolidx] / osimnorm$bwht[r] 
      }
      
    }
    
  }

}
    

# write to file
write_csv(osimnorm, file.path(srcfolder, "force_sdp_results_updated_normalised_sdp.csv"))
