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

# paleobuddy 
library(paleobuddy)

# ape
library(ape)

# ggplot
library(ggplot2)


###
# set parameters for accuracy simulations

# speciation
lambda <- list(0.3, 
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
           c(0.2, 0.05),
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
                         c("DT", "UT", "DT")),
                   rep(30, 7))
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
               "low_unc", "high_unc",
               "m1")
key <- key[, c(5, 1:4)]

# add booleans to make sim code easier
key$model <- c(rep(2, 63), 1)
key$unc <- c(rep("mid", 61), "low", "high", "mid")
key$diff_bins <- c(rep(FALSE, 59), 5, 2, rep(FALSE, 3))

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

# making time bins
make_bins <- function(age, unc) {
  # find range based on unc
  range <- switch(unc,
                  "mid" = {c(age / 6, age / 5)},
                  "low" = {c(age / 10, age / 8)},
                  "high" = {c(age / 4, age / 3)})
  
  # create bins vector
  bins <- c(0)
  
  # while bins doesn't have age
  while (!(age %in% bins)) {
    # if the highest bin is less than the max range, just add age
    if (max(bins) + range[2] > age) {
      bins <- c(bins, age)
    } else {
      # if not, draw a new bin
      bin <- runif(1, range[1], range[2])
      
      # if max(bins) + bin is higher than age, add age instead
      if (max(bins) + bin > age) {
        bins <- c(bins, age)
      } else {
        # add to bins
        bins <- c(bins, max(bins) + bin)
        
      }
    }
  }
  
  # return bins
  return(bins)
}

###
# write auxiliary simulation functions

# simulate one rep
simulate_rep <- function(rates, age, shifts,
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
  pShifts <- shifts[[3]]
  
  # nFinal based on coverage
  nFinal <- ifelse(coverage, 5, 10)
  
  # get bins
  if (coverage) {
    # if coverage, just the shifts (include the final sim time)
    bins <- c(shifts[[1]], age)
  } else {
    # if not, make bins
    bins <- make_bins(age, unc)
  }
  
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
                                             rShifts = pShifts,
                                             bins = bins,
                                             returnAll = TRUE))
    
    # make sure at least some species got sampled
    while (length(unique(fossils$Species)) < 1) {
      fossils <- suppressMessages(sample.clade(sim, psi, tMax, 
                                               rShifts = pShifts,
                                               bins = bins,
                                               returnAll = TRUE))
    }
    
    # start specimens and ranges data frame
    specimens <- data.frame(matrix(nrow = 0, ncol = 4))
    ranges <- data.frame(matrix(nrow = 0, ncol = 3))

    # loop through species
    for (i in 1:length(unique(fossils$Species))) {
      # which species is this
      sp <- unique(fossils$Species)[i]
      
      # get vector of occurrences for that species
      occs <- fossils[fossils$Species == sp, -which(colnames(fossils) == "SampT")]
      
      # and get true occurrence times
      true_occs <- fossils[fossils$Species == sp, c("Species", "Extant", "SampT")]
      
      # check if it is extant
      ext_sp <- fossils$Extant[fossils$Species == sp][1]
      
      # check model
      if (model == 1) {
        # add all occurrences to specimens
        specimens <- rbind(specimens, occs)
      } else {
        # check if there are multiple unique occurrences
        if (nrow(occs) > 1 && length(unique(occs$MaxT)) > 1) {
          # add first and last occurrences to specimens
          specimens <- rbind(specimens, occs[1, ], occs[nrow(occs), ])
        } else {
          # add just one occurrence
          specimens <- rbind(specimens, occs[1, ])
        }
      }
      
      # get true range
      range <- c(max(true_occs$SampT),
                 ifelse(sum(true_occs$Extant) > 0, 0, min(true_occs$SampT)))
      
      # add to ranges
      ranges <- rbind(ranges, c(sp, range))
    }
    
    # reorder and rename columns
    specimens <- specimens[, c(1, 4, 3, 2)]
    colnames(specimens) <- c("taxon", "min_age", "max_age", "status")
    
    # change status column to extant and extinct
    specimens$status <- c("extinct", "extant")[specimens$status + 1]
    
    # name columns for ranges
    colnames(ranges) <- c("taxon", "first_age", "last_age")
    
    # conditions (depending if it's a coverage sim or not)
    if (coverage) {
      # we want to have between 5 and 20 species
      cond <- length(unique(specimens$taxon)) >= 5 && 
        length(unique(specimens$taxon)) <= 50
      
      # if cond is false, redraw everything
      if (!cond) {
        # draw age
        age <- runif(1, 5, 15)
        tMax <- age
        
        # draw rates
        lambda <- rlnorm(3, -2, 0.5)
        mu <- rlnorm(3, -2, 0.5)
        psi <- rlnorm(3, -0.125, 0.5)
        
        # shifts
        lShifts <- mShifts <- pShifts <- seq(0, age, age/3)[-4]
        
        # bins
        bins <- c(lShifts, age)
      }
    }
    else {
      cond <- length(unique(specimens$taxon)) > 10 &&
        length(unique(specimens$taxon)) < 500 &&
        length(unique(specimens$taxon)) / length(sim$TS) < 0.9 &&
        length(unique(specimens$taxon)) / length(sim$TS) > 0.2
    }
  }
  
  # record true values for cov sims
  true_vals <- c(lambda, mu, psi, age)
  
  # return sim, ranges and k
  return(list(SIM = sim, SPECIMENS = specimens, RANGES = ranges, 
              TV = true_vals, BINS = bins))
}

# simulate one set
simulate_set <- function(n_key, reps, rates, age, base_dir,
                         model = 2, unc = FALSE, extant_singletons = FALSE,
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
        shifts[[i]] <- seq(0, age, age / n_stages[i])[-(n_stages[i] + 1)]
      }
    }
    
    # get bins to bin fossil data on
    bins <- seq(0, age, age / 3)
    if (max(n_stages) == 2) bins <- seq(0, age, age / 2)
    
    # if we're on the case where we need different bins
    if (diff_bins) {
      # set those bins
      bins <- seq(0, age, age / diff_bins)
    }
  }
  
  # create vectors for number of species, sampled, and percentage sampled
  n_sp <- c()
  n_sampled <- c()
  perc_sampled <- c()
  
  # create list for sim, fossil ranges, and k
  sims <- vector("list", reps)
  
  # create seeds - reps*100 apart to ensure a lot of possible seeds
  if (coverage) seeds <- runif(reps, 0, 69 * reps * 100) else 
    seeds <- runif(reps, (n_key - 1) * reps * 100, n_key * reps * 100)
  
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
      psi <- rlnorm(3, -0.125, 0.5)
      
      # make rates list
      rates <- list(lambda = lambda,
                    mu = mu,
                    psi = psi)
      
      # shifts list
      shifts <- rep(list(seq(0, age, age / 3)[-4]), 3)
    }
    
    # run sim
    sim_rep <- simulate_rep(rates, age, shifts, 
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
                  col.names = FALSE, row.names = FALSE, 
                  quote = FALSE, sep = "\t")
      
      # add true values to true_vals
      true_vals <- rbind(true_vals, sim_rep$TV)
    }
    
    # get sim, specimens and ranges
    sim <- sim_rep$SIM
    specimens <- sim_rep$SPECIMENS
    ranges <- sim_rep$RANGES

    # calculate numbers of interest
    n_sp <- c(n_sp, length(sim$TS))
    n_sampled <- c(n_sampled, length(unique(specimens$taxon)))
    perc_sampled <- c(perc_sampled, n_sampled[rep]/n_sp[rep])
    
    # save specimens
    smart_dir_create(paste0(base_dir, "/specimens"))
    write.table(specimens, paste0(base_dir, "specimens/taxa_", rep, ".tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    
    # save true ranges
    smart_dir_create(paste0(base_dir, "/true_ranges"))
    write.table(ranges, paste0(base_dir, "true_ranges/ranges_", rep, ".tsv"),
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
    
    # get highest number of rates for this set
    max_len <- max(unlist(lapply(1:length(rates), 
                                 function(x) length(rates[[x]]))))
    
    # get make true rates by expanding rates that are not of len max_len
    true_rates <- lapply(1:length(rates), function(x) {
      if (length(rates[[x]]) < max_len) {
        rep(rates[[x]][1], max_len)
      } else {
        rates[[x]]
      }
    })
    
    # name them
    names(true_rates) <- c("lambda", "mu", "psi")
    
    # make it a data frame
    true_rates <- as.data.frame(true_rates)
    
    # write it to file
    write.table(true_rates, paste0(base_dir, "true_vals.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    # if it is, write true values data frame
    
    # name columns
    colnames(true_vals) <- c("lambda1", "lambda2", "lambda3",
                             "mu1", "mu2", "mu3", 
                             "psi1", "psi2", "psi3", "age")
    
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
                       2, FALSE, FALSE, 0, coverage = TRUE)
  
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

# ##
# # plot time for each sim
# create times vector
# times <- c()
# 
# # iterate through all model 3 analyses of 100 gens
# for (i in 1:63) {
#   # read times vector
#   times_set <- read.table(paste0("/Users/petrucci/Documents/research/",
#                                  "skyfbdr_simstudy/dont_commit/time_",
#                                  i, ".tsv"))
# 
#   # get the total time for this sim set
#   times_sum <- sum(times_set)
# 
#   # add it to times
#   times <- c(times, times_sum)
# }
# 
# # get the total times
# times_total <- 10E6/10 * times / 60 / 60
# 
# # add it to n_sampled
# n_sampled$time <- c(times_total, rep(NA, 4))
# 
# ggplot(n_sampled[-c(64:67), ], aes(y = time)) +
#   geom_histogram() +
#   theme_bw()
