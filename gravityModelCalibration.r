# ================================================================
# Name    :   Calibrate parameters for the lognormal distribution
# Author  :   Rawle Prince
# Date    :   07/12/2015
# Version :   1.2
# =================================================================
require(minpack.lm)
require(MASS)

# -----------------------------------------------------------------
# Implements a version of sumIfs as in excel, with the restriction 
# that the arrays must have the same lengths. There is also an 
# additionanl parameter for dividing the result. This can be useful 
# for normalising the result (i.e. so it represents a proportion) 
# as in the case in which it is currently used
# -----------------------------------------------------------------
adj_sumifs <- function(sum_array, condition_array, f, m=1){
	n <- length(condition_array)
	sm = 0
	if (n == length(condition_array)){
	  fun <- function(x,i){if (f (condition_array[i])){sum_array[i] + x}else{x} }
	  sm <- Reduce(fun,1:n,0)
	} 
	ifelse(m <= 0, sm , sm/m) 
}

# -----------------------------------------------------------------
# estimate the proportions of parameters in between distance bands
# -----------------------------------------------------------------
estimate.inrange <- function(vals,dist,lower,upper,total){
  n <- length(lower)
  if (n == length(upper)){
    sapply(1:n, function(i){ ifelse(i < n ,
                                adj_sumifs(vals,dist, (function(x) x >= lower[i] && x < upper[i]),total) ,
                                adj_sumifs(vals,dist, (function(x) x >= lower[i]) , total)
                             ) }
          )
  }else{
    as.numeric()
  }
}


# -----------------------------------------------------------------
#               REGULARISATION
#
# The procedure works well when the lengths of the low, up and 
# target files are less that or equal to 160. If they are not, the
# performance is very bad. These functions reset these vectors into 
# the range of good performance by extending the interval of low 
# and high, and re-estimating the target values, so that their 
# lenghts are no greater that 160. 
# -----------------------------------------------------------------
intervals.max <- function(vals, mx, maxItvWdt = 10000, lnLim = 2000000){
  ln <- length(vals)
  vs <- vals
  if (ln > mx && ln <= lnLim) {
      for (j in 2 : maxItvWdt ){
        vs <- sapply(seq(from = 1, to = ln, by = j), function(i) vals[i]   )
        if (length(vs) <= mx) break
      }
  } else if (ln > lnLim) {
    simpleWarning(paste("The input vectors are greater than", lnLim,  "and are not reguarised.", sep = " "))
  }
  vs
}

# -----
# perform the regularisation
# -----
regularise.inputs <- function(lower, upper, targets, normTgt, lim = 160){
  ln <- length(lower)
  if (ln == length(upper) && ln == length(targets)){
    if (ln <= lim){
      new_targets <- targets
      if (normTgt){
        sm <- sum(targets)
        new_targets <- sapply(targets, function(x){ if(sm <= 0){x}else{x/ sm }})
      }
      
      list(low = lower, up = upper, targets = new_targets)
    }else{
      new_low <- intervals.max (lower, lim)
      n <- length(new_low)
      new_up  <- c(new_low[2:n], new_low[n] + (new_low[n] - new_low[n-1]) )
      new_targets <- estimate.inrange(vals = targets,dist = lower  ,lower = new_low , upper = new_up,total = sum(targets))
      list(low = new_low, up = new_up, targets = new_targets)
    }
  }else{
    simpleError("The lengths of the input vectors are not the same")
  }
}


# ----
# creates a stratified sample of the OD file based on cost factors.
# The number of stratas is determined by the proportion of the
# cost matrix relative to the maximum process length, then a 
# a random sample of each strata is taken proportationaltely 
# ----
strat_sample_od_matrix <- function(old_od_matrix, od_ln, max_len){
  temp_od_matrix <- old_od_matrix
  df <- data.frame()
  num_stratas <-  round(1.5 + od_ln/max_len)  
  costs <- old_od_matrix$Distance  
  mcost <- round(max(costs))
  temp_od_matrix$Stratas <- cut(costs, seq(0,  mcost , by = mcost / num_stratas ) )
  
  # proportionately sample each strata
  for (cost_factor in unique(temp_od_matrix$Stratas) ){
    if (!is.na(cost_factor) ){
      od_sample <-  subset(temp_od_matrix , Stratas == cost_factor)
      ln <- length(od_sample[,1])
      rn <- length(od_sample[1,])
      od_sample <- od_sample[sample(1 : ln ,  max(1, round(ln * max_len / od_ln))   , replace = FALSE, prob = NULL),  ]
      df <-  rbind(df, od_sample)
    }
  }
  df
}


# -----
#                   FITTING THE LOG NORMAL DISTRIBUTION
# It is assumed that the demand values are always at the last 
# column of the trip length distribution file
# -----
fit.lognormalDist <- function(fileName){
  Data_lognorm<- read.csv(fileName, header = FALSE, strip.white = TRUE, sep = ",")
  m <- length(Data_lognorm[,1])
  lognorm_mtx <- as.matrix(Data_lognorm[2:m,])
  
  low <- as.numeric( lognorm_mtx[,1])
  up <- as.numeric(lognorm_mtx[,2])
  val <- as.numeric(lognorm_mtx[,4])
  vs <- Reduce(function(cs, i) c(cs, as.numeric(replicate( round(val[i]) , (low[i] + up[i]) /2) ))  ,1:(m-1),as.numeric())
  
  Model_lognorm <- fitdistr(vs, 'lognormal')
  estm <- Model_lognorm$estimate
  #write.table(Model_lognorm$estimate, "lognormal_estimates.csv", sep=",")
  out_mtx <- matrix(c(estm[1],estm[2],0,0,0,0), ncol=6)
  
  write.matrix(out_mtx, "lognormal_estimates.csv", sep =",")
}

# ----
# estimate deterrence parameters for lognormal distributions
# ----
calculate_Det_ptns <- function(alpha, beta, pxa, low,up, distances, eF){
  temp <- as.numeric()
  if ( length(pxa) == length(distances) && length(low) == length(up) )
  {
    ln_values <- as.numeric(Map(function(pa,d) eF * pa * dlnorm(d, meanlog = alpha, sdlog = beta),pxa,distances))
    temp <- estimate.inrange (ln_values,distances,low,up, total = sum(ln_values))
  }
  temp
}

# -----
# main function for computing gravity parameters:
# 1 -   Read the files, extract the lower and upper bounds, the target from the target file
# 2 -  lower bounds are at the first column of the target file,
#       upper bounds are at the seond column and target values are
#       at the last column of the target file.
# -----
compute_gravity_parms <- function(source_file, trg_file, expF, mu, sd, 
                                  fn = "logNrm", targets.norm = TRUE, pnt = 0.4){
  
  od_file <-read.csv(source_file, header = TRUE, strip.white = TRUE, sep = ",")
  target_file <- read.csv(trg_file, header = TRUE, strip.white = TRUE, sep = ",")
  
  ln <- length(od_file[,1])
  print(paste("original length: ", ln))
   
  # take a stratified sample of the trip distribution matrix by costs
  # (or dsitance, etc) if the number of inter-zonal costs exceed 220000
  od_file <- strat_sample_od_matrix(od_file, ln , 215000) 

  low <- Filter(function(x) !is.na(x), target_file[,1])
  up  <- Filter(function(x) !is.na(x), target_file[,2])
  m <- length(target_file[1,])
  target <- Filter(function(x) !is.na(x),target_file[,m])
  
  # Sometimes there is an extra column in the target file, in which
  # case target = as.numeric() and the process will fail. In such 
  # instances, we shall asssume that the actual values are at the 
  # proceeding column and select those values as the target
  if (length(target) == 0 && m > 1){
    target <- Filter(function(x) !is.na(x),target_file[,m - 1])
  } 
  
  target <- target[1: length(low)]
  
  # For performance reasons, we shall pre-preocess all vectors before
  # before processing. There are form of pre-processing:
  # --
  # 1. For those target (also upper and lower bounds) vectors of length
  # less that 220, the default option is to normalise the values in the 
  # target vector so that it correspons with proportions. 
  # --
  # 2.If the length of the target vector is greater than 220, the intervals
  # (defined in upper and lower) are increased, and new target proportions
  # are calculated, until the length of the target vector.
  # 
  regularised <- regularise.inputs(low, up, target, normTgt = targets.norm)
  low <- regularised$low 
  up <- regularised$up  
  target <- regularised$targets

  #--------------------------------------------------------------------------------------
  # Calculate Aij x Bij for the purposes, each of which is in adjacent columns. 
  # Note the output from EMME gives: Aij = Oi x Aij  and Bij = Dj x Bij
  #--------------------------------------------------------------------------------------
  dist <- od_file$Distance
  
  od_file <- subset(od_file, select = c(Ai, Bj))
  od_matrix <- as.matrix(od_file)
  
  num_rows  <- length(od_matrix[,1])
  odab <- sapply(1:num_rows, function(i){ od_matrix[i,1] %*% od_matrix[i,2]})
   
  print(paste("length targets: " , length(target)))
  print(paste("length low: " , length(low)))
  print(paste("length up: " , length(up)))
  print(paste("length odab: " , length(odab)))
  print(paste("length distance: " , length(dist)))
  
  calibrated_Prms <- nlsLM(target  ~ calculate_Det_ptns(alpha = a,beta = b, pxa = odab,low = low, up = up, distances = dist, eF = expF),
                        start = list(a = mu, b = sd ), 
                        lower=c(0.00001,0.00001),
                        upper=c(100000.0,100.0),
                        control = nls.lm.control(ftol = 0.001, 
                                                 ptol = 0.001, 
                                                 gtol = 0, diag = list(), epsfcn = 0,
                                                 factor = 100, maxfev = integer(), maxiter = 15, nprint = 0),
                        trace = T)
  
  #mu <- summary(calibrated_Prms)$coefficients[1, 1]
  #std <- summary(calibrated_Prms)$coefficients[2, 1]
  result <- summary(calibrated_Prms)$coefficients
  result_vals <- matrix(c(result[1,1],result[2,1],result[1,2] , result[2,2], result[1,3],result[2,3]), ncol=6)
  write.matrix(result_vals, "lognormal_estimates.csv", sep =",")
  #write.table(c(mu,std) , "lognormal_estimates.csv", sep=",")
}

########################################################################################
## call from batch file with arguments
# 1 argument- the path to the input file
# 6 arguments - the path to the source & target files, together with the
#              starting values for mu, sd, and the expansion factor
########################################################################################
args <- commandArgs(trailingOnly = TRUE)# get some command line arguments
len <- length(args)

if (len == 1){
  # calculating the starting values
  trg_file   <- file.path("../User_Inputs" , as.character(args[1]))
  fit.lognormalDist(trg_file)
}else if(len ==  6){ 
  # calibrating the parameters
  ef <- as.numeric(args[3])
  mu <- as.numeric(args[4])
  sd <- as.numeric(args[5])
  n <- as.numeric(args[6])
  source_file   <- file.path("../Processing" , as.character(args[1]) )
  trg_file   <- file.path("../User_Inputs" , as.character(args[2]))
  compute_gravity_parms(source_file, trg_file, expF = ef, mu, sd, fn = "logNrm")
  
} else {
  print("error calling args. See below:  ")
  print(args)
}