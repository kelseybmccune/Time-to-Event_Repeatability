##############################################
# Combine simulation results (on HPC/katana) #
##############################################

home.wd <- "/srv/scratch/z5394590/survival_repeatability/"

# initialise the result storage
results <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "results/raw", sep = "/"))

# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]

# inner loop through individual sim-by-model files
for (job in res.list) {
  # get the job number
  job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
  # load in the individual simulation results
  load(job)
  # rbind to results data frame
  results <- rbind(results, res)
  # remove the current results to ensure no duplicates
  rm(res)
}


# save the single result data frame within the results folder
setwd(home.wd)
save(list = "results", file = "collated_sim_results.RDATA")



##### Combine all RDATA file of dataframes into one list

home.wd <- "/srv/scratch/z5394590/survival_repeatability/"

# set path to simulated data folder
setwd(paste(home.wd, "results/data", sep = "/"))

# get a list of all Rdata files in the folder
rdata_files <- list.files()[grepl("simdat_", list.files(), fixed = T)]

# initialise empty lists to store the simulated data frames
dat_list <- list()
exploded_dat_q_list <- list()

# loop through each Rdata file and load the data frames and combine into lists
for (file in rdata_files) {
  load(file)
  if (exists("dat")) {
    dat_list[[file]] <- dat
  }
  if (exists("exploded_dat_q")) {
    exploded_dat_q_list[[file]] <- exploded_dat_q
  }
}

# Save the lists of data frames into a single Rdata file
save(dat_list, exploded_dat_q_list, file = "collated_sim_data.RDATA")
