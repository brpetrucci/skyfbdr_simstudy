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

# ggplot
library(ggplot2)

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
key_extra <- cbind(rbind(c("DT", "C", "C"),
                         c("DT", "UT", "C"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT"),
                         c("DT", "UT", "DT")),
                   rep(30, 10))
key_extra <- as.data.frame(key_extra)
colnames(key_extra) <- colnames(key)
key <- rbind(key, as.data.frame(key_extra))

# name each sim set
key$names <- c(paste0(rep(c("base", "lLI3", "lLI2", 
                            "lLD3", "lLD2", "lUT", "lDT",
               "mLI3", "mLI2", "mLD3", "mLD2", "mUT", "mDT",
               "pLI3", "pLI2", "pLD3", "pLD2", "pUT", "pDT")), 
               "_", c(rep(10, 19), rep(20, 19), rep(30, 19))),
               "free", "not_free",
               "stages_5", "stages_2",
               "no_extant_singles", "yes_extant_singles",
               "m1_no_unc", "m1_unc", "m2_no_unc", "m2_unc")
key <- key[, c(5, 1:4)]

# add booleans to make sim code easier
key$model <- c(rep(3, 63), 1, 1, 2, 2)
key$unc <- c(rep(FALSE, 64), TRUE, FALSE, TRUE)
key$extant_singles <- c(rep(FALSE, 62), TRUE, rep(FALSE, 4))
key$diff_bins <- c(rep(FALSE, 59), 5, 2, rep(FALSE, 6))

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
                         model, unc, extant_singletons,
                         coverage = FALSE) {
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
  
  # nFinal based on coverage
  nFinal <- ifelse(coverage, 5, 10)
  
  ##
  # run simulations
  
  # create conditions boolean
  cond <- FALSE
  
  # run until we fulfill conditions
  while (!cond) {
    # run BD simulation - make sure we get 10+ species
    sim <- bd.sim(n0, lambda, mu, tMax, 
                  lShifts = lShifts, mShifts = mShifts, nFinal = c(nFinal, Inf))
    
    # run fossil sampling
    fossils <- suppressMessages(sample.clade(sim, psi, tMax, 
                                             rShifts = pShifts))
    
    # make sure at least some species got sampled
    while (length(unique(fossils$Species)) < 1) {
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
      colnames(ranges) <- c("taxon", "max_age", "min_age")
    }
    
    # get the number of species sampled in each interval
    k_sampled <- colSums(k_nums)
    
    # conditions (depending if it's a coverage sim or not)
    if (coverage) {
      cond <- nrow(ranges) >= 5 && nrow(ranges) <= 20
      
      # if cond is false, redraw everything
      if (!cond) {
        # draw age
        age <- runif(1, 5, 15)
        
        # draw rates
        lambda <- rlnorm(3, -2, 0.5)
        mu <- rlnorm(3, -2, 0.5)
        psi <- rexp(3, 1)
        
        # shifts
        lShifts <- mShifts <- pShifts <- seq(0, age, age/3)[-4]
        
        # bins
        bins <- c(lShifts, age)
      }
    }
    else {
      cond <- nrow(ranges) > 10 &&
        nrow(ranges) < 500 &&
        nrow(ranges)/length(sim$TS) < 0.9 &&
        nrow(ranges)/length(sim$TS) > 0.2 &&
        all(k_sampled > 0)
    }
  }
  
  # record true values for cov sims
  true_vals <- c(lambda, mu, psi, age)
  
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
  
  # return sim, ranges and k
  return(list(SIM = sim, RANGES = ranges, K = k, TV = true_vals, BINS = bins))
}

# simulate one set
simulate_set <- function(n_key, reps, rates, age, base_dir,
                         model = 3, unc = FALSE, extant_singletons = FALSE,
                         diff_bins = 0, coverage = FALSE) {
  # if it is not a coverage simulation set, can prepare rates as normal
  if (!coverage) {
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
    
    # if we're on the case where we need different bins
    if (diff_bins) {
      # set those bins
      bins <- seq(0, age, age/diff_bins)
    }
  }
  
  # create vectors for number of species, sampled, and percentage sampled
  n_sp <- c()
  n_sampled <- c()
  perc_sampled <- c()
  
  # create list for sim, fossil ranges, and k
  sims <- vector("list", reps)
  
  # create seeds - reps*100 apart to ensure a lot of possible seeds
  if (coverage) seeds <- runif(reps, 0, 69*reps*100) else 
    seeds <- runif(reps, (n_key - 1)*reps*100, n_key*reps*100)
  
  # save seeds
  save(seeds, file = paste0(base_dir, "seeds.RData"))
  
  # start true values data frame if coverage
  true_vals <- data.frame(matrix(nrow = 0, ncol = 4))
  
  # iterate through reps
  for (rep in 1:reps) {
    # print some info
    print(paste0("key: ", n_key, " rep: ", rep, 
                 " seed: ", seeds[rep]))
    
    # set seed 
    set.seed(seeds[rep])
    
    # if it is a coverage sim set, gotta draw the values
    if (coverage) {
      # draw age
      age <- runif(1, 5, 15)
      
      # draw rates
      lambda <- rlnorm(3, -2, 0.5)
      mu <- rlnorm(3, -2, 0.5)
      psi <- rexp(3, 1)
      
      # make rates list
      rates <- list(lambda = lambda,
                    mu = mu,
                    psi = psi)
      
      # shifts list
      shifts <- rep(list(seq(0, age, age/3)[-4]), 3)
      
      # bins
      bins <- c(shifts[[1]], age)
    }
    
    # run sim
    sim_rep <- simulate_rep(rates, age, shifts, bins, 
                            model = model, unc = unc, 
                            extant_singletons = extant_singletons,
                            coverage = coverage)
    
    # if coverage, need to same timeline and add to true_vals
    if (coverage) {
      # get final bins
      bins <- sim_rep$BINS
      
      # write timeline
      smart_dir_create(paste0(base_dir, "/times"))
      write.table(t(bins[-c(1, length(bins))]), 
                  paste0(base_dir, "/times/times_", rep, ".tsv"),
                  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      
      # add true values to true_vals
      true_vals <- rbind(true_vals, sim_rep$TV)
    }
    
    # get sim and samples
    sim <- sim_rep$SIM
    ranges <- sim_rep$RANGES
    k <- sim_rep$K
    
    # calculate numbers of interest
    n_sp <- c(n_sp, length(sim$TS))
    if (extant_singletons) {
      ranges_fs <- ranges[ranges$max_age != 0 | ranges$min_age != 0, ]
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
  
  # if not a coverage sim
  if (!coverage) {
    # write timeline
    write.table(t(bins[-c(1, length(bins))]), 
                paste0(base_dir, "times.tsv"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    # if it is, write true values data frame
    
    # name columns
    colnames(true_vals) <- c("lambda", "mu", "psi", "age")
    
    # save true_vals
    write.table(true_vals, 
                paste0(base_dir, "/true_vals.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # make numbers into a data frame
  nums <- data.frame(n_sp = n_sp, n_sampled = n_sampled,
                     perc_sampled = perc_sampled)
  
  # return numbers of interest
  return(nums)
}

# write final accuracy simulate function
simulate_acc <- function(key, reps, lambda, mu, psi, reps_dir) {
  # create refs directory
  refs_dir <- paste0(reps_dir, "refs/")
  smart_dir_create(refs_dir)
  
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
    diff_bins_set <- key_set$diff_bins
    
    # create base directory
    base_dir <- paste0(reps_dir, key_set$names, "/")
    smart_dir_create(base_dir)
    
    # run simulation set
    nums <- simulate_set(n_key, reps, rates_set, age_set, base_dir,
                         model_set, unc_set, extant_singletons_set,
                         diff_bins_set, coverage = FALSE)

    # add to data frames
    min_nums <- rbind(min_nums, colMins(nums))
    max_nums <- rbind(max_nums, colMaxes(nums))
    mean_nums <- rbind(mean_nums, colMeans(nums))
    
    # start lines to add to refs script
    refs_lines <- vector("character", 8)
    
    # add directory name
    refs_lines[1] <- paste0('dir <- "', key_set$names, '"')
    
    # check number of rates to run analysis with
    n_lambda <- ifelse(key_set$lambda == "C", 1,
                       ifelse(key_set$lambda %in% c("LI2", "LD2"), 2, 3))
    n_mu <- ifelse(key_set$mu == "C", 1,
                   ifelse(key_set$mu %in% c("LI2", "LD2"), 2, 3))
    n_psi <- ifelse(key_set$psi == "C", 1,
                    ifelse(key_set$psi %in% c("LI2", "LD2"), 2, 3))
    
    # correct it if we're on the free or not_free sims
    if (key_set$names == "free") {
      n_mu <- n_psi <- 3
    }
    if (key_set$names == "not_free") {
      n_mu <- 1
    }
    
    # correct it if we're on the stages_5 or stages_2 sims
    if (key_set$names == "stages_5") {
      n_lambda <- n_mu <- n_psi <- 5
    }
    if (key_set$names == "stages_2") {
      n_lambda <- n_mu <- n_psi <- 2
    }
    
    # add lines for lambda, mu and psi number
    refs_lines[2] <- paste0("n_lambda <- ", n_lambda)
    refs_lines[3] <- paste0("n_mu <- ", n_mu)
    refs_lines[4] <- paste0("n_psi <- ", n_psi)
    
    # add lines for age and other parameters
    refs_lines[5] <- paste0("age <- ", age_set)
    refs_lines[6] <- paste0("model <- ", model_set)
    refs_lines[7] <- paste0("unc <- ", unc_set)
    refs_lines[8] <- paste0("ext_singles <- ", extant_singletons_set)
    
    # write lines to refs directory
    writeLines(refs_lines, paste0(refs_dir, "ref_", n_key, ".Rev"))
  }
  
  # save key to reps_dir
  write.table(key, paste0(reps_dir, "key.tsv"), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # name nums data frames
  colnames(min_nums) <- colnames(max_nums) <- colnames(mean_nums) <-
    c("n_sp", "n_sampled", "perc_sampled")
  
  # return the data frames
  return(list(MINS = min_nums, MAXES = max_nums, MEANS = mean_nums))
}

# create function for simulating coverage
simulate_cov <- function(reps, reps_dir) {
  # run simulation set
  nums <- simulate_set("cov", reps, NULL, NULL, reps_dir,
                       3, FALSE, FALSE, 0, coverage = TRUE)
  
  # make nums a data frame
  nums <- as.data.frame(nums)
  
  # name nums data frames
  colnames(nums) <- c("n_sp", "n_sampled", "perc_sampled")
  
  # return the data frames
  return(nums)
}

###
# run simulations - coverage

# create reps directory
cov_reps_dir <- paste0("/Users/petrucci/Documents/research/skyfbdr_simstudy/",
                       "simulation/coverage/")
smart_dir_create(cov_reps_dir)

# set number of reps
reps <- 1000

# run sims
nums_cov <- simulate_cov(reps, cov_reps_dir)

# save nums_cov
write.table(nums_cov, paste0("cov_reps_dir", "nums_cov.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

###
# run simulations - accuracy

# create reps directory
acc_reps_dir <- paste0("/Users/petrucci/Documents/research/skyfbdr_simstudy/",
                       "simulation/accuracy/")
smart_dir_create(acc_reps_dir)

# save key
write.table(key, 
            paste0(acc_reps_dir, "key.tsv"),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# set number of reps
reps <- 100

# run sims
nums_acc <- simulate_acc(key, reps, lambda, mu, psi, acc_reps_dir)

# save nums_acc
write.table(nums_acc$MEANS, paste0("acc_reps_dir", "means_acc.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(nums_acc$MAXES, paste0("acc_reps_dir", "maxes_acc.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(nums_acc$MINS, paste0("acc_reps_dir", "mins_acc.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# get n_sampled data frames
n_sampled <- data.frame(mean = nums_acc$MEANS$n_sampled, 
                        min = nums_acc$MINS$n_sampled,
                        max = nums_acc$MAXES$n_sampled,
                        age = key$age,
                        name = key$names,
                        lambda = key$lambda)
n_sampled_10 <- n_sampled[n_sampled$age == 10, ]
n_sampled_20 <- n_sampled[n_sampled$age == 20, ]
n_sampled_30 <- n_sampled[n_sampled$age == 30, ]

# plot each n_sampled
ggplot(n_sampled_10, aes(x = 1:nrow(n_sampled_10), y = mean,
                         col = lambda)) +
  geom_point() +
  geom_segment(aes(x = 1:nrow(n_sampled_10), y = min, yend = max)) +
  theme_bw()
ggplot(n_sampled_20, aes(x = 1:nrow(n_sampled_20), y = mean,
                         col = lambda)) +
  geom_point() +
  geom_segment(aes(x = 1:nrow(n_sampled_20), y = min, yend = max)) +
  theme_bw()
ggplot(n_sampled_30, aes(x = 1:nrow(n_sampled_30), y = mean,
                         col = lambda)) +
  geom_point() +
  geom_segment(aes(x = 1:nrow(n_sampled_30), y = min, yend = max)) +
  theme_bw()
