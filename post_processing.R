# for ess
library(coda)

# for plotting
library(reshape2)
library(ggplot2)

# base directory to save on (reps will be names taxa_rep.tsv)
base_dir <- "/Users/petrucci/Documents/research/skyfbdr_simstudy/"

# function to read a certain quantiles and check if true value is within
check_credible <- function(log, trueVal, quants) {
  # get quantiles
  unlist(lapply(1:length(quants), function(x) {
    qs <- quantile(log, c(0.5 - quants[x]/2, 0.5 + quants[x]/2))
    trueVal > qs[1] & trueVal < qs[2]
  }))
}

# function to get back ess and credible intervals for variables
get_vals <- function(vars, log, trueVals, quants) {
  # create dataframes to hold values
  creds <- data.frame(matrix(nrow = length(vars), ncol = length(quants)))
  errors <- c()
  ess <- c()
  
  # iterate through vars
  for (i in 1:length(vars)) {
    # get true value
    var_true <- trueVals[vars[i]]
    
    # add to creds
    creds[i, ] <- check_credible(log[, vars[i]], var_true, quants)
    
    # and to errors and ess
    errors <- c(errors, abs(mean(log[, vars[i]]) - var_true) / var_true)
    ess <- c(ess, effectiveSize(log[, vars[i]]))
  }
  
  # name stuff
  colnames(creds) <- quants
  rownames(creds) <- names(errors) <- names(ess) <- vars
  
  # return the data frames
  return(list(CREDS = creds, ERRORS = unlist(errors), ESS = ess))
}

# create function to evaluate each model's output
validate <- function(sim_kind) {
  # number of reps
  nreps <- 1000
  
  # true values
  trueVals <- read.delim(paste0(base_dir, "simulation/", 
                                sim_kind, "/true_vals.tsv"), sep = "\t")
  
  # quantiles of interest
  quants <- seq(0.05, 0.95, 0.05)
  
  # create data frames for counts of true/false
  lambda1_creds <- lambda2_creds <- lambda3_creds <-
    mu1_creds <- mu2_creds <- mu3_creds <- 
    psi1_creds <- psi2_creds <- psi3_creds <-
    data.frame(matrix(nrow = 0, ncol = length(quants)))
  errors <- data.frame(matrix(nrow = 0, ncol = 9))
  ess <- data.frame(matrix(nrow = 0, ncol = 9))
  
  # read log files
  for (rep in 1:nreps) {
    # get log
    log <- read.delim(paste0(base_dir, "output/", sim_kind, "/", 
                             "skyfbdr_", 
                             rep, ".log"))[, 5:13]
    
    # reorder log
    log <- log[, c(3, 2, 1, 6, 5, 4, 9, 8, 7)]
    
    # just to make the names a bit easier
    colnames(log) <- sort(gsub("\\.", "", colnames(log)))
    
    # burnin - 25%
    log <- log[(nrow(log)/10):nrow(log), ]
    
    # get creds, errors, and ess for each variable
    vals <- get_vals(colnames(log), log, trueVals[rep, ], quants)
    creds <- vals$CREDS
    errors_rep <- vals$ERRORS
    ess_rep <- vals$ESS
    
    # add to creds data frames
    lambda1_creds <- rbind(lambda1_creds, creds[1, ])
    lambda2_creds <- rbind(lambda2_creds, creds[2, ])
    lambda3_creds <- rbind(lambda3_creds, creds[3, ])
    mu1_creds <- rbind(mu1_creds, creds[4, ])
    mu2_creds <- rbind(mu2_creds, creds[5, ])
    mu3_creds <- rbind(mu3_creds, creds[6, ])
    psi1_creds <- rbind(psi1_creds, creds[7, ])
    psi2_creds <- rbind(psi2_creds, creds[8, ])
    psi3_creds <- rbind(psi3_creds, creds[9, ])
    
    errors <- rbind(errors, errors_rep)
    ess <- rbind(ess, ess_rep)
  }
  
  # make creds just one dataframe
  creds_full <- t(data.frame(lambda1 = colSums(lambda1_creds) / nreps,
                             lambda2 = colSums(lambda2_creds) / nreps,
                             lambda3 = colSums(lambda3_creds) / nreps,
                             mu1 = colSums(mu1_creds) / nreps,
                             mu2 = colSums(mu2_creds) / nreps,
                             mu3 = colSums(mu3_creds) / nreps,
                             psi1 = colSums(psi1_creds) / nreps,
                             psi2 = colSums(psi2_creds) / nreps,
                             psi3 = colSums(psi3_creds) / nreps))
  
  # make rownames into column
  creds_full_melted <- melt(creds_full)
  colnames(creds_full_melted) <- c("rate", "cred_int", "prop_within")
  
  # check min ESS
  print(paste0("Mininum ESS: ", min(ess)))
  
  colnames(errors) <- colnames(ess) <- c("lambda1", "lambda2", "lambda3",
                                         "mu1", "mu2", "mu3",
                                         "psi1", "psi2", "psi3")
  
  # plot
  plot <- ggplot(creds_full_melted, aes(x = cred_int, 
                                        y = prop_within, color = rate)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "Credible interval width", y = "Truth recovery proportion",
         color = "Rate") +
    theme_bw()
  
  print(plot)
  
  return(list(ESS = ess, ERR = errors, CREDS = creds_full))
}

cov_val <- validate("coverage")
ess <- cov_val$ESS
err <- cov_val$ERR
creds <- cov_val$CREDS
