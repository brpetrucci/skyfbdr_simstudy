######################################################
##                    Chapter 3                     ##
##      Treeless skyline FBDR simulation study      ##
##                Simulation script                 ##
##            Bruno do Rosario Petrucci             ##
######################################################

###
# load packages

# devtools
library(devtools)

# paleobuddy - loading since I want to include local changes
load_all()

# ape
library(ape)

###
# set parameters for accuracy simulations

# speciation
lambda<- list(0.3, 
              c(0.1, 0.2, 0.4),
              c(0.1, 0.4),
              c(0.4, 0.2, 0.1),
              c(0.4, 0.1),
              c(0.1, 0.4, 0.1),
              c(0.4, 0.1, 0.4))

# extinction
mu <- list(0.1,
           c(0.05, 0.1, 0.2),
           c(0.05, 0.2),
           c(0.2, 0.1, 0.05), 
           c(0.05, 0.2),
           c(0.05, 0.2, 0.05),
           c(0.2, 0.05, 0.2))

# fossil sampling
psi <- lambda

# name them
names(lambda) <- names(mu) <- names(psi) <-
  c("C", "LI3", "LI2", "LD3", "LD2", "UT", "DT")

# set key to define which values we want to match
key <- data.frame(lambda = c(names(lambda),
                             rep("C", 12)),
                  mu = c(rep("C", 7), names(mu)[-1],
                         rep("C", 6)),
                  psi = c(rep("C", 13), names(psi)[-1]))

# repeat it 3 times to add origin ages
key <- cbind(rbind(key, key, key), 
             age = c(rep(10, 19), rep(20, 19), rep(30, 19)))

# add extra sims
key_extra <- cbind(rbind(c("C", "LD2", "DT"),
                         c("LI2", "UT", "C"),
                         c("DT", "C", "LI2"),
                         c("DT", "C", "C"),
                         c("DT", "UT", "C"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT")),
                   rep(30, 11))
key_extra <- as.data.frame(key_extra)
colnames(key_extra) <- colnames(key)
key <- rbind(key, as.data.frame(key_extra))

# name each sim set
key$names <- c(paste0(rep(c("base", "lLI3", "lLI2", 
                            "lLD3", "lLD2", "lUT", "lDT",
               "mLI3", "mLI2", "mLD3", "mLD2", "mUT", "mDT",
               "pLI3", "pLI2", "pLD3", "pLD2", "pUT", "pDT")), 
               "_", c(rep(10, 19), rep(20, 19), rep(30, 19))),
               "mult_stages123", "mult_stages231", "mult_stages312",
               "free", "not_free",
               "no_extant_singles", "yes_extant_singles",
               "m1_no_unc", "m1_unc", "m2_no_unc", "m2_unc")
key <- key[, c(5, 1:4)]

# add booleans to make sim code easier
key$model <- c(rep(3, 64), 1, 1, 2, 2)
key$unc <- c(rep(FALSE, 65), TRUE, FALSE, TRUE)
key$extant_singles <- c(rep(FALSE, 64), TRUE, rep(FALSE, 3))

###
# write minor auxiliary functions

# smarter dir.create function
smart_dir_create <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# equivalent to colSums for max
colMaxes <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, max(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

# same for min
colMins <- function(df) {
  # make return
  res <- c()
  
  # iterate through columns
  for (i in 1:ncol(df)) {
    # add max to res
    res <- c(res, min(df[, i]))
  }
  
  # name res
  names(res) <- colnames(df)
  
  # return res
  return(res)
}

###
# write auxiliary simulation functions

# simulate one rep
simulate_rep <- function(rates, age, shifts, bins,
                         model, unc, extant_singletons) {
  # need to run sim, run sampling, then make sampling into ranges
  # k should be 0 and 1 based on bins, or just count for each bin
  # if extant singletons is true, include unsampled extant species
  
  ## set parameters
  # number of initial species
  n0 <- 1
  
  # total simulation time
  tMax <- age
  
  # lambda
  lambda <- rates[[1]]
  
  # mu
  mu <- rates[[2]]
  
  # psi
  psi <- rates[[3]]
  
  # shifts
  lShifts <- shifts[[1]]
  mShifts <- shifts[[2]]
  pShifts <- shifts[[2]]
  
  ##
  # run simulations
  
  # create conditions boolean
  cond <- FALSE
  
  # run until we fulfill conditions
  while (!cond) {
    # run BD simulation - make sure we get 10+ species
    sim <- bd.sim(n0, lambda, mu, tMax, 
                  lShifts = lShifts, mShifts = mShifts, nFinal = c(10, Inf))
    
    # run fossil sampling
    fossils <- suppressMessages(sample.clade(sim, psi, tMax, 
                                             rShifts = pShifts))
    
    # make sure at least some species got sampled
    while (length(unique(fossils$Species)) < 5) {
      fossils <- suppressMessages(sample.clade(sim, psi, tMax, 
                                               rShifts = pShifts))
    }
    
    # start ranges data frame
    ranges <- data.frame(matrix(nrow = 0,
                                ncol = ifelse(unc, 5, 3)))
    # if uncertainty is true, need a column for min and max for FAD and LAD
    
    # fossil counts
    k <- data.frame(matrix(nrow = 0, ncol = (length(bins) - 1)))

    # loop through species
    for (i in 1:length(unique(fossils$Species))) {
      # which species is this
      sp <- unique(fossils$Species)[i]
      
      # get vector of occurrences for that species
      occs <- fossils$SampT[fossils$Species == sp]
      
      # find fad and lad
      fad <- max(occs)
      lad <- min(occs)
      
      # start a vector for k for this species
      k_sp <- c()
      
      # iterate through bins
      for (j in 1:(length(bins) - 1)) {
        # add the number of occurrences in this bin to k_sp
        k_sp <- c(k_sp, sum(occs > bins[j] & occs <= bins[j + 1]))
      }
      
      # if there is uncertainty, bin FAD and LAD
      if (unc) {
        fad_max <- min(bins[bins >= fad])
        fad_min <- max(bins[bins < fad])
        
        lad_max <- min(bins[bins >= lad])
        lad_min <- max(bins[bins < lad])
        
        # add to ranges
        ranges <- rbind(ranges, c(sp, fad_max, fad_min, lad_max, lad_min))
      }
      else {
        # just add fad and lad to ranges
        ranges <- rbind(ranges, c(sp, fad, lad))
      }
      
      # if model is 3, make k 0s and 1s
      if (model == 3) {
        k_sp <- ifelse(k_sp == 0, 0, 1)
      }
      
      # add k_sp to k
      k <- rbind(k, k_sp)
      
      # save just the number part of k
      k_nums <- k
    }
    
    # if model is 1 or 2, collapse k into a vector of interval fossil counts
    if (model != 3) {
      k <- c("t", colSums(k))
      
      # name k
      names(k) <- c("taxon", paste0("int_", 1:(length(bins) - 1)))
    }
    else {
      # add taxon names
      k <- cbind(unique(fossils$Species), k)
      
      # name k
      colnames(k) <- c("taxon", paste0("int_", 1:(length(bins) - 1)))
    }
    
    # if uncertainty is true, name things differently
    if (unc) {
      colnames(ranges) <- c("taxon", "fad_max", "fad_min", "lad_max", "lad_min")
    }
    else {
      colnames(ranges) <- c("taxon", "fad", "lad")
    }
    
    # if extant_singletons is true, add extant singletons, if any, to ranges
    if (extant_singletons) {
      # find which taxa are extant
      ext <- paste0("t", which(sim$EXTANT))
      
      # find which are singletons
      ext_singles <- ext[!(ext %in% fossils$Species)]
      
      # add them to ranges with 0 for every age, if there are any
      if (length(ext_singles) > 0) {
        # make a zeros data frame
        zeros_df <- data.frame(matrix(0, nrow = length(ext_singles),
                                      ncol = ncol(ranges) - 1))
        
        # add taxa names to it
        ext_singles_df <- cbind(ext_singles, zeros_df)
        
        # name it
        colnames(ext_singles_df) <- colnames(ranges)
        
        # add new columns to range
        ranges <- rbind(ranges, ext_singles_df)
      }
    }
    
    # get the number of species sampled in each interval
    k_sampled <- colSums(k_nums)
    
    # check conditions
    cond <- length(unique(fossils$Species)) > 10 &&
              length(unique(fossils$Species)) < 500 &&
                length(unique(fossils$Species))/length(sim$TS) < 0.9 &&
                  length(unique(fossils$Species))/length(sim$TS) > 0.2 &&
                    all(k_sampled > 0)
  }
  
  # return sim, ranges and k
  return(list(SIM = sim, RANGES = ranges, K = k))
}

# simulate one set
simulate_set <- function(n_key, reps, rates, age, base_dir,
                         model = 3, unc = FALSE, extant_singletons = FALSE) {
  # create shifts list
  shifts <- vector("list", 3)
  
  # and number of stages per rate
  n_stages <- c()
  
  # iterate through rates to see how many shifts per rate
  for (i in 1:length(rates)) {
    # check number of stages for this rate
    n_stages <- c(n_stages, length(rates[[i]]))
    
    # make shifts vector
    if (n_stages[i] > 1) { 
      shifts[[i]] <- seq(0, age, age/n_stages[i])[-(n_stages[i] + 1)]
    }
  }
  
  # get bins to bin fossil data on
  bins <- seq(0, age, age/3)
  if (max(n_stages) == 2) bins <- seq(0, age, age/2)
  
  # create vectors for number of species, sampled, and percentage sampled
  n_sp <- c()
  n_sampled <- c()
  perc_sampled <- c()
  
  # create list for sim, fossil ranges, and k
  sims <- vector("list", reps)
  
  # create seeds - reps*100 apart to ensure a lot of possible seeds
  seeds <- runif(reps, (n_key - 1)*reps*100, n_key*reps*100)
  
  # save seeds
  save(seeds, file = paste0(base_dir, "seeds.RData"))
  
  # iterate through reps
  for (rep in 1:reps) {
    # print some info
    print(paste0("key: ", n_key, " rep: ", rep, 
                 " seed: ", seeds[rep]))
    
    # set seed 
    set.seed(seeds[rep])
    
    # run sim
    sim_rep <- simulate_rep(rates, age, shifts, bins, 
                            model = model, unc = unc, 
                            extant_singletons = extant_singletons)
    
    # get sim and samples
    sim <- sim_rep$SIM
    ranges <- sim_rep$RANGES
    k <- sim_rep$K
    
    # calculate numbers of interest
    n_sp <- c(n_sp, length(sim$TS))
    if (extant_singletons) {
      ranges_fs <- ranges[ranges$fad != 0 | ranges$lad != 0, ]
      n_sampled <- c(n_sampled, length(unique(ranges_fs$taxon)))
    } else {
      n_sampled <- c(n_sampled, length(unique(ranges$taxon)))
    }
    perc_sampled <- c(perc_sampled, n_sampled[rep]/n_sp[rep])
    
    # save ranges
    smart_dir_create(paste0(base_dir, "/ranges"))
    write.table(ranges, paste0(base_dir, "ranges/taxa_", rep, ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # save k
    smart_dir_create(paste0(base_dir, "/fossil_counts"))
    if (!is.data.frame(k)) k <- t(k)
    write.table(k, paste0(base_dir, "fossil_counts/k_", rep, ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # append to list of simulations
    sims[[rep]] <- sim
  }
  
  # save sim list
  save(sims, file = paste0(base_dir, "sim_list.RData"))
  
  # write timeline
  write.table(t(bins[-c(1, length(bins))]), 
              paste0(base_dir, "timeline", rep, ".tsv"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # make numbers into a data frame
  nums <- data.frame(n_sp = n_sp, n_sampled = n_sampled,
                     perc_sampled = perc_sampled)
  
  # return numbers of interest
  return(nums)
}

# write final simulate function
simulate_all <- function(key, reps, lambda, mu, psi, reps_dir) {
  # create data frame to hold numbers of interest
  min_nums <- max_nums <- mean_nums <- data.frame(matrix(nrow = 0, ncol = 3))
  
  # iterate through key
  for (n_key in 1:nrow(key)) {
    # get key values 
    key_set <- key[n_key, ]
    
    # make rates 
    rates_set <- list(lambda = lambda[[key_set$lambda]], 
                      mu = mu[[key_set$mu]],
                      psi = psi[[key_set$psi]])
    
    # get age 
    age_set <- as.numeric(key_set$age)
    
    # get other values
    model_set <- as.numeric(key_set$model)
    unc_set <- key_set$unc
    extant_singletons_set <- key_set$extant_singles
    
    # create base directory
    base_dir <- paste0(reps_dir, key_set$names, "/")
    smart_dir_create(base_dir)
    
    # run simulation set
    nums <- simulate_set(n_key, reps, rates_set, age_set, base_dir,
                         model_set, unc_set, extant_singletons_set)

    # add to data frames
    min_nums <- rbind(min_nums, colMins(nums))
    max_nums <- rbind(max_nums, colMaxes(nums))
    mean_nums <- rbind(mean_nums, colMeans(nums))
  }
  
  # name nums data frames
  colnames(min_nums) <- colnames(max_nums) <- colnames(mean_nums) <-
    c("n_sp", "n_sampled", "perc_sampled")
  
  # return the data frames
  return(list(MINS = min_nums, MAXES = max_nums, MEANS = mean_nums))
}

###
# run simulations

# create reps directory
reps_dir <- paste0("/Users/petrucci/Documents/research/skyfbdr_simstudy/",
                   "simulation/accuracy/replicates/")
smart_dir_create(reps_dir)

# set number of reps
reps <- 100

# run sims
nums <- simulate_all(key, reps, lambda, mu, psi, reps_dir)
