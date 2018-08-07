# Copyright (C) 2018 Eduard Grebe, Alex Welte and Stellenbosch University
#
# Some code re-used from inctools: <https://github.com/SACEMA/inctools>.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

P <- function(t, parameters) {
  1 / (1 + exp( -(parameters[1] + parameters[2]*t + parameters[3]*t^2 + 
                    parameters[4]*t^3) ) )
}

dPda <- function(t, parameters) {
  ( exp(parameters[1] + parameters[2]*t + parameters[3]*t^2 + parameters[4]*t^3) * 
      (parameters[2] + 2*parameters[3]*t + 3* parameters[4]*t^2) ) / 
    ( 1 + exp(parameters[1] + parameters[2]*t + parameters[3]*t^2 + parameters[4]*t^3)  )^2
}

dPdt <- function(P, lambda = 0) {
  return(-P * lambda)
}

PR <- function(t, parameters) {
  exp(-exp(parameters[1] + (parameters[2]) * log(t)))
}

IR <- function(prevH, prevR, mdri, frr, bigT) {
  prevH * (prevR - frr)/((1 - prevH) * (mdri - frr * bigT))
}

sample_frac_groups <- function(tbl, size, replace = FALSE, weight=NULL) {
  import::from(magrittr, "%>%")
  # regroup when done
  grps <- tbl %>% dplyr::groups() %>% base::unlist() %>% base::as.character()
  # check length of groups non-zero
  keep <- tbl %>% dplyr::summarise() %>% dplyr::sample_frac(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% dplyr::right_join(keep, by=grps) %>% dplyr::group_by_(grps)
}

bootstrap_params <- function(data, sex = "Female", age.range = c(15,35), 
                             n_bootstraps = 12500, cores = 4, gc_step = 1000,
                             debug = FALSE) {
  if(debug) {browser()}
  
  import::from(magrittr, "%>%")
  import::from(foreach, "%dopar%")
  
  if (sex == "both") {
    data_fit <- dplyr::filter(data, 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else if (sex == "Female") {
    data_fit <- dplyr::filter(data, sex == "Female", 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else if (sex == "Male") {
    data_fit <- dplyr::filter(data, sex == "Male", 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else {stop("Sex must be 'both', 'Female' or 'Male'")}
  #browser()
  wards <- unique(data_fit$ward)
  
  ward_data_grouped <- list()
  for (i in 1:length(wards)) {
    ward_data_grouped[[i]] <- data_fit %>%
      dplyr::filter(ward == wards[i]) %>%
      dplyr::group_by(cluster)
  }
  
  # options for glm2
  tolerance <- 1e-08
  maxit <- 50000
  
  cluster <- parallel::makeCluster(cores, type = "FORK", outfile = "")
  doParallel::registerDoParallel(cluster)
  pb <- utils::txtProgressBar(min = 1, max = n_bootstraps, style = 3)
  library(foreach)
  gc_iterations <- seq(1,n_bootstraps,gc_step)
  gc()
  params <- foreach::foreach(i = 1:n_bootstraps, .combine = rbind, .inorder = FALSE, 
                             .packages = "dplyr",
                             .errorhandling = "remove") %dopar% 
                             {
                               if (i %in% gc_iterations) {gc()}
                               if(!exists("pb")) pb <- utils::txtProgressBar(min = 1, max = n, style = 3)
                               
                               ward_data_resampled <- lapply(ward_data_grouped,
                                                             function(df){
                                                               df %>%
                                                                 sample_frac_groups(1, replace = TRUE) %>%
                                                                 dplyr::ungroup()
                                                             })
                               bootstrapped_data <- dplyr::bind_rows(ward_data_resampled)
                               
                               fit_prev <- glm2::glm2(formula = hiv ~ 1 + 
                                                        I(age_years) +
                                                        I(age_years^2) + 
                                                        I(age_years^3), 
                                                      family = stats::binomial(link = "logit"),
                                                      data = bootstrapped_data, 
                                                      control = stats::glm.control(epsilon = tolerance, 
                                                                                   maxit = maxit,
                                                                                   trace = FALSE))
                               
                               fit_prevR <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                                                         I(log(age_years)),
                                                       family = stats::binomial(link = "cloglog"), 
                                                       data = bootstrapped_data,
                                                       control = stats::glm.control(epsilon = tolerance,
                                                                                    maxit = maxit, 
                                                                                    trace = FALSE))
                               
                               paramsvec <- c(fit_prev$coefficients, fit_prevR$coefficients)
                               utils::setTxtProgressBar(pb, i)
                               return(paramsvec)
                             }
  close(pb)
  parallel::stopCluster(cluster)
  return(params)
}

extract_incidences <- function(params, deriv_t = 0, mort.spline, 
                               age.range = c(15,35), age.step = 1, cores = 4, 
                               gc_step = 1000, mdri, frr, bigT, debug = FALSE) {
  if (debug) {browser()}
  import::from(foreach, "%dopar%")
  ages <- seq(age.range[1], age.range[2], age.step)
  
  cluster <- parallel::makeCluster(cores, type = "FORK", outfile = "")
  doParallel::registerDoParallel(cluster)
  pb <- utils::txtProgressBar(min = 1, max = nrow(params), style = 3)
  n <- nrow(params)
  library(foreach)
  gc_iterations <- seq(1,n,gc_step)
  gc()
  incarray <- foreach::foreach(i = 1:n, 
                               .combine = function(...) abind::abind(..., along=3), 
                               .inorder = FALSE, 
                               .export = c("P","dPda","dPdt","PR","IR"),
                               .packages = "abind",
                               .errorhandling = "stop") %dopar% 
                               {
                                 if (i %in% gc_iterations) {gc()}
                                 if(!exists("pb")) pb <- utils::txtProgressBar(min = 1, max = n, style = 3)
                                 iterate_mat <- matrix(nrow = length(ages), ncol = 4)

                                 iterate_mat[,1] <- P(t = ages, parameters = unname(as.matrix(params[i,1:4])))
                                 
                                 iterate_mat[,2] <- (1/(1-iterate_mat[,1])) * (dPda(t = ages, parameters = unname(as.matrix(params[i,1:4]))) + dPdt(a = ages, deriv = deriv_t)) + # (1/(1-p(a))) * (dPda + dPdt)
                                   iterate_mat[,1] * predict(mort.spline, x = ages)$y # p(a)* mort(a)
                                 
                                 
                                 iterate_mat[,3] <- PR(t = ages, parameters = unname(as.matrix(params[i,5:6])))
                                 iterate_mat[,4] <- IR(prevH = iterate_mat[,1], 
                                                       prevR = iterate_mat[,3], 
                                                       mdri = mdri, 
                                                       frr = frr, 
                                                       bigT = bigT)
                                 # new version end
                                 rownames(iterate_mat) <- ages
                                 colnames(iterate_mat) <- c("PrevH","IncPrev","PrevR","IncR")
                                 utils::setTxtProgressBar(pb, i)
                                 return(iterate_mat)
                               }
  close(pb)
  parallel::stopCluster(cluster)
  dimnames(incarray)[[3]] <- paste0("iteration.",1:dim(incarray)[[3]])
  return(incarray)
}

incidence <- function(data, sex = "Female", age.range = c(15,35), mort.spline, 
                      incmat,  deriv_t = 0, mdri, rse_mdri, frr, rse_frr, bigT, 
                      alpha = 0.05, debug = FALSE) {
  
  if(debug) {browser()}
  
  # options for glm2
  tolerance <- 1e-08
  maxit <- 50000
  
  # Function for sigmaW
  sigmaW <- function(W,sigma1,sigma2,rho) {
    sqrt(W^2*sigma1^2 + (1-W)^2*sigma2^2 + 2*rho*sigma1*sigma2*W*(1-W))
  }
  
  if (sex == "both") {
    data_fit <- dplyr::filter(data, 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else if (sex == "Female") {
    data_fit <- dplyr::filter(data, sex == "Female", 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else if (sex == "Male") {
    data_fit <- dplyr::filter(data, sex == "Male", 
                              age_years >= get("age.range")[1], 
                              age_years <= get("age.range")[2])
  } else {stop("Sex must be 'both', 'Female' or 'Male'")}
  
  
  fit_prev <- glm2::glm2(formula = hiv ~ 1 + 
                           I(age_years) +
                           I(age_years^2) + 
                           I(age_years^3), 
                         family = stats::binomial(link = "logit"),
                         data = data_fit, 
                         control = stats::glm.control(epsilon = tolerance, 
                                                      maxit = maxit,
                                                      trace = FALSE))
  fit_prevR <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                            I(log(age_years)),
                          family = stats::binomial(link = "cloglog"), 
                          data = data_fit,
                          control = stats::glm.control(epsilon = tolerance,
                                                       maxit = maxit, 
                                                       trace = FALSE))
  paramsvec <- unname(c(fit_prev$coefficients, fit_prevR$coefficients))
  
  
  IncWmat <- matrix(nrow = nrow(incmat[,,1]), ncol = 15)
  colnames(IncWmat) <- c("Age","PrevH","PrevH.sigma","IncPrev","IncPrev.sigma","PrevR","PrevR.sigma","IncR","IncR.sigma","Rho","W","IncW","IncW.sigma","IncW.LB","IncW.UB")
  IncWmat[,1] <- as.numeric(as.character(rownames(incmat[,,1])))
  n_estimates <- dim(incmat)[3]
  
  IncWmat[,2] <- P(t = IncWmat[,1], parameters = paramsvec[1:4]) # PrevH
  # PrevH.sigma -> loop
  
  # MAJOR BUG: 1/1-p normalisation by number of susceptibles was left out!!
  IncWmat[,4] <- (1/(1-IncWmat[,2])) * ( dPda(t = IncWmat[,1], parameters = paramsvec[1:4]) + dPdt(a = IncWmat[,1], deriv = deriv_t)) + IncWmat[,2] * predict(mort.spline, x = IncWmat[,1])$y # IncPrev
  # IncPrev.sigma -> loop
  IncWmat[,6] <- PR(t = IncWmat[,1], parameters = paramsvec[5:6]) # PrevR
  # PrevR.sigma - loop
  IncWmat[,8] <- IR(prevH = IncWmat[,2], 
                    prevR = IncWmat[,6], 
                    mdri = mdri, 
                    frr = frr, 
                    bigT = bigT) # IncR
  # IncR.sigma -> loop
  # Rho, W, IncW, IncW.sigma, IncW.LB, IncW.UB -> loop
  
  for (i in 1:nrow(IncWmat)) {
    IncWmat[i,3] <- sd(incmat[i,1,]) # PrevH.sigma
    IncWmat[i,7] <- sd(incmat[i,3,]) # PrevR.sigma
    
    sigma1 <- sd(incmat[i,2,]) # convenience var for sigma(IncPrev)
    IncWmat[i,5] <- sigma1 #IncPrev.sigma
    
    inc_temp <- inctools::incprops(PrevH = IncWmat[i,2],
                                   RSE_PrevH = 0,
                                   PrevR = IncWmat[i,6],
                                   RSE_PrevR = 0,
                                   Boot = FALSE,
                                   MDRI = mdri*365.25,
                                   RSE_MDRI = rse_mdri,
                                   FRR = frr,
                                   RSE_FRR = rse_frr,
                                   BigT = bigT*365.25,
                                   Covar_HR = 0) 
    sigma.MDRI.FRR <- as.numeric(as.character(inc_temp$Incidence.Statistics$RSE)) * 
      as.numeric(as.character(inc_temp$Incidence.Statistics$Incidence)) # for sigma on IncR from MDRI and FRR
    sigma2 <- sqrt(sd(incmat[i,4,])^2 + sigma.MDRI.FRR^2) # total sigma for IncR
    IncWmat[i,9] <- sigma2 #IncR.sigma
    
    # Find optimal weighting factor
    rho <- cor(x = incmat[i,1,], y = incmat[i,3,])
    IncWmat[i,10] <- rho
    
    W <- ( sigma2^2 - rho*sigma1*sigma2 ) / (sigma1^2 + sigma2^2 - 2*rho*sigma1*sigma2)
    if (W>1) {W <- 1}
    if (W<0) {W <- 0}
    IncWmat[i,11] <- W
    
    # IncPrev = IncWmat[i,4] # convenience variabke
    # IncR = IncWmat[i,8] # convenience variable
    IncWmat[i,12] <- W*IncWmat[i,4] + (1-W)*IncWmat[i,8] # IncW (total incidence)
  }
  IncWmat[,13] <- sigmaW(W = IncWmat[,11], sigma1 = IncWmat[,5], sigma2 = IncWmat[,9], rho = IncWmat[,10]) # IncW.sigma
  IncWmat[,14] <- IncWmat[,12] - qnorm(1-(alpha/2))*IncWmat[,13] # IncW.LB
  IncWmat[,14] <- ifelse(IncWmat[,14]<0,0,IncWmat[,14]) # set lower bound to zero if negative
  IncWmat[,15] <- IncWmat[,12] + qnorm(1-(alpha/2))*IncWmat[,13] # IncW.UB
  IncWmat[,15] <- ifelse(IncWmat[,15]<0,0,IncWmat[,15]) # set upper bound to zero if negative
  IncWmat[,12] <- ifelse(IncWmat[,12]<0,0,IncWmat[,12]) # set point estimate to zero if negative
  
  return(tibble::as_tibble(IncWmat))
}

average_incidence <- function(data, sex = "Female", age.range.fit = c(15,35), step = 0.25, age.range.av = c(15,35), 
                              deriv_t = 0, 
                              weighting = "sampling", poptable = NULL, mort.spline.f, mort.spline.m, mdri, 
                              frr, bigT, Wvec.m, Wvec.f, cores = 4, n_bootstraps = 10000, alpha = 0.05, debug = FALSE)
{
  if (debug) {browser()}
  
  # Internal functions
  int <- function(func,ll,ul,...) {
    return(
      cubature::adaptIntegrate(f = func, 
                               lowerLimit = ll, 
                               upperLimit = ul, 
                               ...)$integral
    )
  }
  
  func_mult <- function(x,func1,func2) {
    return(func1(x) * func2(x))
  }
  
  incidence_weighted <- function(incfunc,densityfunc,ll,ul) {
    integral <- int(func = func_mult, ll = ll, ul = ul, func1 = incfunc, func2 = densityfunc)
    normalisation <- int(func = densityfunc, ll = ll, ul = ul)
    incidence <- integral/normalisation
    return(incidence)
  }
  
  incidence_weighted_bysex <- function(incfunc.m, incfunc.f, densityfunc.m, densityfunc.f,ll,ul) {
    Wm <- int(densityfunc.m,ll,ul)
    Wf <- int(densityfunc.f,ll,ul)
    Im <- incidence_weighted(incfunc = incfunc.m, densityfunc = densityfunc.m, ll, ul)
    If <- incidence_weighted(incfunc = incfunc.f, densityfunc = densityfunc.f, ll, ul)
    incidence <- ( (Wm * Im) + (Wf * If) ) / (Wm + Wf)
  }
  
  # This function's code
  if (weighting == "population" & is.null(poptable)) {stop("If weighting by population, a population table must be provided")}
  import::from(magrittr, "%>%")
  import::from(foreach, "%dopar%")
  
  # options for glm2
  tolerance <- 1e-08
  maxit <- 50000
  
  age_min <- get("age.range.fit")[1]
  age_max <- get("age.range.fit")[2]
  age_vec <- seq(age_min,age_max,step)
  ll <- age.range.av[1]
  ul <- age.range.av[2]
  
  data_fit.f <- dplyr::filter(data, sex == "Female", 
                              age_years >= age_min, 
                              age_years <= age_max)
  data_fit.m <- dplyr::filter(data, sex == "Male", 
                              age_years >= age_min, 
                              age_years <= age_max)
  
  
  # Point estimate
  
  if (sex != "Male") {
    fit_prev.f <- glm2::glm2(formula = hiv ~ 1 + 
                               I(age_years) +
                               I(age_years^2) + 
                               I(age_years^3), 
                             family = stats::binomial(link = "logit"),
                             data = data_fit.f, 
                             control = stats::glm.control(epsilon = tolerance, 
                                                          maxit = maxit,
                                                          trace = FALSE))
    
    fit_prevR.f <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                                I(log(age_years)),
                              family = stats::binomial(link = "cloglog"), 
                              data = data_fit.f,
                              control = stats::glm.control(epsilon = tolerance,
                                                           maxit = maxit, 
                                                           trace = FALSE))
    
    PrevH.f <- P(t = age_vec, parameters = unname(fit_prev.f$coefficients))
    
    # MAJOR BUG: 1/1-p normalisation by number of susceptibles was left out!!
    IncPrev.f <- (1 / (1 - PrevH.f)) * ( dPda(t = age_vec, parameters = unname(fit_prev.f$coefficients))  + dPdt(a = age_vec, deriv = deriv_t)) + 
      PrevH.f * predict(mort.spline.f, x = age_vec)$y
    PrevR.f <- PR(t = age_vec, parameters = unname(fit_prevR.f$coefficients))
    IncR.f  <- IR(prevH = PrevH.f, 
                  prevR = PrevR.f, 
                  mdri = mdri, 
                  frr = frr, 
                  bigT = bigT)
    IncW.f <- Wvec.f*IncPrev.f + (1-Wvec.f)*IncR.f
    inc.spline.f <- splinefun(x = age_vec, y = IncW.f)
  }
  
  if (sex != "Female") {
    fit_prev.m <- glm2::glm2(formula = hiv ~ 1 + 
                               I(age_years) +
                               I(age_years^2) + 
                               I(age_years^3), 
                             family = stats::binomial(link = "logit"),
                             data = data_fit.m, 
                             control = stats::glm.control(epsilon = tolerance, 
                                                          maxit = maxit,
                                                          trace = FALSE))
    
    fit_prevR.m <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                                I(log(age_years)),
                              family = stats::binomial(link = "cloglog"), 
                              data = data_fit.m,
                              control = stats::glm.control(epsilon = tolerance,
                                                           maxit = maxit, 
                                                           trace = FALSE))
    
    PrevH.m <- P(t = age_vec, parameters = unname(fit_prev.m$coefficients))
    
    # MAJOR BUG: 1/1-p normalisation by number of susceptibles was left out!!
    IncPrev.m <- (1 / (1 - PrevH.m)) * ( dPda(t = age_vec, parameters = unname(fit_prev.m$coefficients)) + dPdt(a = age_vec, deriv = deriv_t)) + 
      PrevH.m * predict(mort.spline.m, x = age_vec)$y
    PrevR.m <- PR(t = age_vec, parameters = unname(fit_prevR.m$coefficients))
    IncR.m  <- IR(prevH = PrevH.m, 
                  prevR = PrevR.m, 
                  mdri = mdri, 
                  frr = frr, 
                  bigT = bigT)
    IncW.m <- Wvec.m*IncPrev.m + (1-Wvec.m)*IncR.m
    inc.spline.m <- splinefun(x = age_vec, y = IncW.m)
  }
  
  # Weighting
  # we now create the density weighting functions here, so we can restrict to HIV-negative population after fitting a prevalence model
  if (weighting == "population") {
    
    if (sex != "Male") {
      poptable_fit.f <- poptable %>%
        dplyr::filter(Age >= age_min, Age <= age_max) %>%
        dplyr::select(Age, n = Female) %>%
        dplyr::mutate(n = as.numeric(n))
      
      pop_prev_vec <- P(t = poptable_fit.f$Age, parameters = unname(fit_prev.f$coefficients)) # find prevalences at each age
      poptable_fit.f$n <- poptable_fit.f$n * (1 - pop_prev_vec) # set n to susceptible population
      
      density_spline.f <- splinefun(x = poptable_fit.f$Age, y = poptable_fit.f$n)
    }
    
    if (sex != "Female") {
      poptable_fit.m <- poptable %>%
        dplyr::filter(Age >= age_min, Age <= age_max) %>%
        dplyr::select(Age, n = Male) %>%
        dplyr::mutate(n = as.numeric(n))
      
      pop_prev_vec <- P(t = poptable_fit.m$Age, parameters = unname(fit_prev.m$coefficients)) # find prevalences at each age
      poptable_fit.m$n <- poptable_fit.m$n * (1 - pop_prev_vec) # set n to susceptible population
      
      density_spline.m <- splinefun(x = poptable_fit.m$Age, y = poptable_fit.m$n)
    }
    
  } else if (weighting == "sampling") {
    sampling_density.f <- density(data_fit.f$age_years, na.rm = TRUE)
    sampling_density.m <- density(data_fit.m$age_years, na.rm = TRUE)
    density_spline.f <- splinefun(x = sampling_density.f$x, y = sampling_density.f$y)
    density_spline.m <- splinefun(x = sampling_density.m$x, y = sampling_density.m$y)
  }
  
  #Calculate the average incidence by weighting the incidence functions
  if (sex == "Female") {
    AvInc.PE <- incidence_weighted(incfunc = inc.spline.f, densityfunc = density_spline.f, ll = ll, ul = ul)
  } else if (sex == "Male") {
    AvInc.PE <- incidence_weighted(incfunc = inc.spline.m, densityfunc = density_spline.m, ll = ll, ul = ul)
  } else if (sex == "both") {
    AvInc.PE <- incidence_weighted_bysex(incfunc.m = inc.spline.m, densityfunc.m = density_spline.m, 
                                         incfunc.f = inc.spline.f, densityfunc.f = density_spline.f,
                                         ll = ll, ul = ul)
  } else {stop("Sex must be 'Male', 'Female' or 'both'")}
  
  
  # Bootstrapping for CI
  
  wards <- unique(data$ward)
  ward_data_grouped <- list()
  for (i in 1:length(wards)) {
    ward_data_grouped[[i]] <- data %>%
      dplyr::filter(ward == wards[i]) %>%
      dplyr::group_by(cluster)
  }
  
  cluster <- parallel::makeCluster(cores, type = "FORK", outfile = "")
  doParallel::registerDoParallel(cluster)
  pb <- utils::txtProgressBar(min = 1, max = n_bootstraps, style = 3)
  
  average_incidences <- foreach::foreach(i = 1:n_bootstraps, 
                                         .combine = c, 
                                         .inorder = FALSE,
                                         .errorhandling = "stop") %dopar% 
                                         {
                                           if(!exists("pb")) pb <- utils::txtProgressBar(min = 1, max = n, style = 3)
                                           ward_data_resampled <- lapply(ward_data_grouped,
                                                                         function(df){
                                                                           df %>%
                                                                             sample_frac_groups(1, replace = TRUE) %>%
                                                                             dplyr::ungroup()
                                                                         })
                                           bootstrapped_data <- dplyr::bind_rows(ward_data_resampled) %>%
                                             dplyr::filter(age_years >= age_min, 
                                                           age_years <= age_max)
                                           
                                           bootstrapped_data.f <- bootstrapped_data %>%
                                             filter(sex == "Female")
                                           bootstrapped_data.m <- bootstrapped_data %>%
                                             filter(sex == "Male")
                                           
                                           if (weighting == "sampling") {
                                             sampling_density.f <- density(bootstrapped_data.f$age_years, na.rm = TRUE)
                                             sampling_density.m <- density(bootstrapped_data.m$age_years, na.rm = TRUE)
                                             density_spline.f <- splinefun(x = sampling_density.f$x, y = sampling_density.f$y)
                                             density_spline.m <- splinefun(x = sampling_density.m$x, y = sampling_density.m$y)
                                           }
                                           
                                           if (sex != "Male") {
                                             fit_prev.f <- glm2::glm2(formula = hiv ~ 1 + 
                                                                        I(age_years) +
                                                                        I(age_years^2) + 
                                                                        I(age_years^3), 
                                                                      family = stats::binomial(link = "logit"),
                                                                      data = bootstrapped_data.f, 
                                                                      control = stats::glm.control(epsilon = tolerance, 
                                                                                                   maxit = maxit,
                                                                                                   trace = FALSE),
                                                                      start = unname(fit_prev.f$coefficients))
                                             
                                             fit_prevR.f <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                                                                         I(log(age_years)),
                                                                       family = stats::binomial(link = "cloglog"), 
                                                                       data = bootstrapped_data.f,
                                                                       control = stats::glm.control(epsilon = tolerance,
                                                                                                    maxit = maxit, 
                                                                                                    trace = FALSE),
                                                                       start = unname(fit_prevR.f$coefficients))
                                             
                                             PrevH.f <- P(t = age_vec, parameters = unname(fit_prev.f$coefficients))
                                             # MAJOR BUG: 1/1-p normalisation by number of susceptibles was left out!!
                                             IncPrev.f <- (1 / (1 - PrevH.f)) * ( dPda(t = age_vec, parameters = unname(fit_prev.f$coefficients))  + dPdt(a = age_vec, deriv = deriv_t)) + 
                                               PrevH.f * predict(mort.spline.f, x = age_vec)$y
                                             PrevR.f <- PR(t = age_vec, parameters = unname(fit_prevR.f$coefficients))
                                             IncR.f  <- IR(prevH = PrevH.f, 
                                                           prevR = PrevR.f, 
                                                           mdri = mdri, 
                                                           frr = frr, 
                                                           bigT = bigT)
                                             IncW.f <- Wvec.f*IncPrev.f + (1-Wvec.f)*IncR.f
                                             inc.spline.f <- splinefun(x = age_vec, y = IncW.f)
                                           }
                                           
                                           if (sex != "Female") {
                                             fit_prev.m <- glm2::glm2(formula = hiv ~ 1 + 
                                                                        I(age_years) +
                                                                        I(age_years^2) + 
                                                                        I(age_years^3), 
                                                                      family = stats::binomial(link = "logit"),
                                                                      data = bootstrapped_data.m, 
                                                                      control = stats::glm.control(epsilon = tolerance, 
                                                                                                   maxit = maxit,
                                                                                                   trace = FALSE),
                                                                      start = unname(fit_prev.m$coefficients))
                                             
                                             fit_prevR.m <- glm2::glm2(formula = (1 - recent) ~ 1 + 
                                                                         I(log(age_years)),
                                                                       family = stats::binomial(link = "cloglog"), 
                                                                       data = bootstrapped_data.m,
                                                                       control = stats::glm.control(epsilon = tolerance,
                                                                                                    maxit = maxit, 
                                                                                                    trace = FALSE),
                                                                       start = unname(fit_prevR.m$coefficients))
                                             
                                             PrevH.m <- P(t = age_vec, parameters = unname(fit_prev.m$coefficients))
                                             # MAJOR BUG: 1/1-p normalisation by number of susceptibles was left out!!
                                             IncPrev.m <- (1 / (1 - PrevH.m)) * ( dPda(t = age_vec, parameters = unname(fit_prev.m$coefficients)) + dPdt(a = age_vec, deriv = deriv_t)) + 
                                               PrevH.m * predict(mort.spline.m, x = age_vec)$y
                                             PrevR.m <- PR(t = age_vec, parameters = unname(fit_prevR.m$coefficients))
                                             IncR.m  <- IR(prevH = PrevH.m, 
                                                           prevR = PrevR.m, 
                                                           mdri = mdri, 
                                                           frr = frr, 
                                                           bigT = bigT)
                                             IncW.m <- Wvec.m*IncPrev.m + (1-Wvec.m)*IncR.m
                                             inc.spline.m <- splinefun(x = age_vec, y = IncW.m)
                                           }
                                           
                                           if (sex == "Female") {
                                             average_incidence <- incidence_weighted(incfunc = inc.spline.f, densityfunc = density_spline.f, ll = ll, ul = ul)
                                           } else if (sex == "Male") {
                                             average_incidence <- incidence_weighted(incfunc = inc.spline.m, densityfunc = density_spline.m, ll = ll, ul = ul)
                                           } else if (sex == "both") {
                                             average_incidence <- incidence_weighted_bysex(incfunc.m = inc.spline.m, densityfunc.m = density_spline.m, 
                                                                                           incfunc.f = inc.spline.f, densityfunc.f = density_spline.f,
                                                                                           ll = ll, ul = ul)
                                             
                                           }
                                           utils::setTxtProgressBar(pb, i)
                                           return(average_incidence)
                                         }
  close(pb)
  parallel::stopCluster(cluster)
  
  AvInc.CI <- stats::quantile(average_incidences, probs = c(alpha/2, 1 - alpha/2))
  return(c(PE = AvInc.PE, CI = AvInc.CI))
}
