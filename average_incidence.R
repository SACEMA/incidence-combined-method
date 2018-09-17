# Copyright (C) 2018 Eduard Grebe and Stellenbosch University
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

  setwd("~/dev/incidence-kzn/output/")
  library(feather); library(haven); library(dplyr); library(inctools);
  library(glm2); library(parallel); library(doParallel); library(foreach);
  library(readxl); library(stringr); library(lubridate); library(tidyr);
  library(survey); library(import); library(abind);

# FUNCTIONS

# P <- function(t, parameters) {
#   exp(parameters[1] + parameters[2] * t + parameters[3] * t^2 + parameters[4] *
#         t^3)
# }
#
# dPda <- function(t, parameters) {
#   P(t, parameters) * (parameters[2] + 2*parameters[3]*t + 3*parameters[4]*t^2)
# }
#
# # hard-coded for now
# dPdt <- function(a, deriv = 0) {
#   return(rep(deriv, length(a)))
# }

# New version (logit link)
P <- function(t, parameters) {
  1 / (1 + exp( -(parameters[1] + parameters[2]*t + parameters[3]*t^2 +
                    parameters[4]*t^3) ) )
}

# derivative of logit version
dPda <- function(t, parameters) {
  ( exp(parameters[1] + parameters[2]*t + parameters[3]*t^2 + parameters[4]*t^3) *
      (parameters[2] + 2*parameters[3]*t + 3* parameters[4]*t^2) ) /
    ( 1 + exp(parameters[1] + parameters[2]*t + parameters[3]*t^2 + parameters[4]*t^3)  )^2
}


# hard-coded for now
dPdt <- function(a, deriv = 0) {
  return(rep(deriv, length(a)))
}

PR <- function(t, parameters) {
  exp(-exp(parameters[1] + (parameters[2]) * log(t)))
}

# SIMPLE FORMULA FOR INCIDENCE
IR <- function(prevH, prevR, mdri, frr, bigT) {
  prevH * (prevR - frr)/((1 - prevH) * (mdri - frr * bigT))
}



sample_frac_groups = function(tbl, size, replace = FALSE, weight=NULL) {
  # regroup when done
  grps <- tbl %>% dplyr::groups() %>% base::unlist() %>% base::as.character()
  # check length of groups non-zero
  keep <- tbl %>% dplyr::summarise() %>% dplyr::sample_frac(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% dplyr::right_join(keep, by=grps) %>% dplyr::group_by_(grps)
}

kzn <- read_feather("kzn.feather")
mortality_table <- read_feather("mortality_table.feather")

mortality_f_2013 <- mortality_table %>%
  dplyr::filter(Sex == "Female") %>%
  dplyr::select(age = Age, mort = `2013`)
mort.spline.f <- smooth.spline(x = mortality_f_2013$age, y = mortality_f_2013$mort)

mortality_m_2013 <- mortality_table %>%
  dplyr::filter(Sex == "Male") %>%
  dplyr::select(age = Age, mort = `2013`)
mort.spline.m <- smooth.spline(x = mortality_m_2013$age, y = mortality_m_2013$mort)

mortality_b_2013 <- mortality_table %>%
  dplyr::filter(Sex == "Combined") %>%
  dplyr::select(age = Age, mort = `2013`)
mort.spline.b <- smooth.spline(x = mortality_b_2013$age, y = mortality_b_2013$mort)

poptable <- read_feather("poptable.feather")

load("bs_params_w1535.RData")
load("incarray_w1535.RData")
inc_w1535 <- read_feather("inc_w1535.feather")
women_summary <- read_feather("women_summary.feather")
load("bs_params_m1535.RData")
load("incarray_m1535.RData")
inc_m1535 <- read_feather("inc_m1535.feather")
men_summary <- read_feather("men_summary.feather")

load("bs_params_b1535.RData")
load("incarray_b1535.RData")
inc_b1535 <- read_feather("inc_b1535.feather")
both_summary <- read_feather("both_summary.feather")


# Parameters

mdri <- 224.5231 - 7.71 # Adjustment for NAT screening with DT 100 c/ml, pools of 5 (effective DT 500c/ml) # Subtype C!!
mdriYr <- mdri / 365.25
rse_mdri <- 0.0619132331026123 # NB: uncertainty in NAT adjustment not known or taken into account
frrContext <- 0.00166780137766038
rse_frrContext <- 0.465416630819573

n_bootstraps <- 10000
cores <- 4

averages.samp <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                            Age_L = c(15,20,25,30,15,25,15,15),
                            Age_U = c(20,25,30,35,25,35,30,35),
                            Women = rep(NA,8),
                            Women.LB = rep(NA,8),
                            Women.UB = rep(NA,8),
                            Men = rep(NA,8),
                            Men.LB = rep(NA,8),
                            Men.UB = rep(NA,8),
                            Both = rep(NA,8),
                            Both.LB = rep(NA,8),
                            Both.UB = rep(NA,8))
for (i in 1:nrow(averages.samp)) {
  Women <- average_incidence(data = kzn,
                             sex = "Female",
                             age.range.fit = c(15,35),
                             step = 0.25,
                             age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                             weighting = "sampling",
                             poptable = poptable,
                             mort.spline.f = mort.spline.f,
                             mort.spline.m = mort.spline.m,
                             mdri = mdriYr,
                             frr = frrContext,
                             bigT = 2,
                             Wvec.m = inc_m1535$W,
                             Wvec.f = inc_w1535$W,
                             cores = cores,
                             n_bootstraps = n_bootstraps,
                             alpha = 0.05)
  Men <- average_incidence(data = kzn,
                           sex = "Male",
                           age.range.fit = c(15,35),
                           step = 0.25,
                           age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                           weighting = "sampling",
                           poptable = poptable,
                           mort.spline.f = mort.spline.f,
                           mort.spline.m = mort.spline.m,
                           mdri = mdriYr,
                           frr = frrContext,
                           bigT = 2,
                           Wvec.m = inc_m1535$W,
                           Wvec.f = inc_w1535$W,
                           cores = cores,
                           n_bootstraps = n_bootstraps,
                           alpha = 0.05)
  Both <- average_incidence(data = kzn,
                            sex = "both",
                            age.range.fit = c(15,35),
                            step = 0.25,
                            age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                            weighting = "sampling",
                            poptable = poptable,
                            mort.spline.f = mort.spline.f,
                            mort.spline.m = mort.spline.m,
                            mdri = mdriYr,
                            frr = frrContext,
                            bigT = 2,
                            Wvec.m = inc_m1535$W,
                            Wvec.f = inc_w1535$W,
                            cores = cores,
                            n_bootstraps = n_bootstraps,
                            alpha = 0.05)
  averages.samp$Women[i] <- Women[[1]]
  averages.samp$Women.LB[i] <- Women[[2]]
  averages.samp$Women.UB[i] <- Women[[3]]
  averages.samp$Men[i] <- Men[[1]]
  averages.samp$Men.LB[i] <- Men[[2]]
  averages.samp$Men.UB[i] <- Men[[3]]
  averages.samp$Both[i] <- Both[[1]]
  averages.samp$Both.LB[i] <- Both[[2]]
  averages.samp$Both.UB[i] <- Both[[3]]
}
write_feather(averages.samp, "averages.sampweighted.feather")

averages.pop <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                           Age_L = c(15,20,25,30,15,25,15,15),
                           Age_U = c(20,25,30,35,25,35,30,35),
                           Women = rep(NA,8),
                           Women.LB = rep(NA,8),
                           Women.UB = rep(NA,8),
                           Men = rep(NA,8),
                           Men.LB = rep(NA,8),
                           Men.UB = rep(NA,8),
                           Both = rep(NA,8),
                           Both.LB = rep(NA,8),
                           Both.UB = rep(NA,8))
for (i in 1:nrow(averages.pop)) {
  Women <- average_incidence(data = kzn,
                             sex = "Female",
                             age.range.fit = c(15,35),
                             step = 0.25,
                             age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                             weighting = "population",
                             poptable = poptable,
                             mort.spline.f = mort.spline.f,
                             mort.spline.m = mort.spline.m,
                             mdri = mdriYr,
                             frr = frrContext,
                             bigT = 2,
                             Wvec.m = inc_m1535$W,
                             Wvec.f = inc_w1535$W,
                             cores = cores,
                             n_bootstraps = n_bootstraps,
                             alpha = 0.05)
  Men <- average_incidence(data = kzn,
                           sex = "Male",
                           age.range.fit = c(15,35),
                           step = 0.25,
                           age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                           weighting = "population",
                           poptable = poptable,
                           mort.spline.f = mort.spline.f,
                           mort.spline.m = mort.spline.m,
                           mdri = mdriYr,
                           frr = frrContext,
                           bigT = 2,
                           Wvec.m = inc_m1535$W,
                           Wvec.f = inc_w1535$W,
                           cores = cores,
                           n_bootstraps = n_bootstraps,
                           alpha = 0.05)
  Both <- average_incidence(data = kzn,
                            sex = "both",
                            age.range.fit = c(15,35),
                            step = 0.25,
                            age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                            weighting = "population",
                            poptable = poptable,
                            mort.spline.f = mort.spline.f,
                            mort.spline.m = mort.spline.m,
                            mdri = mdriYr,
                            frr = frrContext,
                            bigT = 2,
                            Wvec.m = inc_m1535$W,
                            Wvec.f = inc_w1535$W,
                            cores = cores,
                            n_bootstraps = n_bootstraps,
                            alpha = 0.05)
  averages.pop$Women[i] <- Women[[1]]
  averages.pop$Women.LB[i] <- Women[[2]]
  averages.pop$Women.UB[i] <- Women[[3]]
  averages.pop$Men[i] <- Men[[1]]
  averages.pop$Men.LB[i] <- Men[[2]]
  averages.pop$Men.UB[i] <- Men[[3]]
  averages.pop$Both[i] <- Both[[1]]
  averages.pop$Both.LB[i] <- Both[[2]]
  averages.pop$Both.UB[i] <- Both[[3]]
}
write_feather(averages.pop, "averages.popweighted.feather")

# For discussion - comparison with Blaizot et al.
average_incidence(data = kzn,
                  sex = "Female",
                  age.range.fit = c(15,35),
                  step = 0.25,
                  age.range.av = c(15,25),
                  weighting = "population",
                  poptable = poptable,
                  mort.spline.f = mort.spline.f,
                  mort.spline.m = mort.spline.m,
                  mdri = mdriYr,
                  frr = frrContext,
                  bigT = 2,
                  Wvec.m = inc_m1535$W,
                  Wvec.f = inc_w1535$W,
                  cores = cores,
                  n_bootstraps = n_bootstraps,
                  alpha = 0.05)
average_incidence(data = kzn,
                  sex = "Female",
                  age.range.fit = c(15,35),
                  step = 0.25,
                  age.range.av = c(25,35),
                  weighting = "population",
                  poptable = poptable,
                  mort.spline.f = mort.spline.f,
                  mort.spline.m = mort.spline.m,
                  mdri = mdriYr,
                  frr = frrContext,
                  bigT = 2,
                  Wvec.m = inc_m1535$W,
                  Wvec.f = inc_w1535$W,
                  cores = cores,
                  n_bootstraps = n_bootstraps,
                  alpha = 0.05)


## *** Synthetic cohort method alone ***

average_incidence_sc <- function(data, sex = "Female", age.range.fit = c(15,35), step = 0.25, age.range.av = c(15,35),
                                 deriv_t = 0,
                                 weighting = "sampling", poptable = NULL, mort.spline.f, mort.spline.m,
                                 cores = 4, n_bootstraps = 10000, alpha = 0.05, debug = FALSE)
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

    PrevH.f <- P(t = age_vec, parameters = unname(fit_prev.f$coefficients))
    IncPrev.f <- (1 / (1 - PrevH.f)) * ( dPda(t = age_vec, parameters = unname(fit_prev.f$coefficients))  + dPdt(a = age_vec, deriv = deriv_t)) +
      PrevH.f * predict(mort.spline.f, x = age_vec)$y
    inc.spline.f <- splinefun(x = age_vec, y = IncPrev.f)
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

    PrevH.m <- P(t = age_vec, parameters = unname(fit_prev.m$coefficients))
    IncPrev.m <- (1 / (1 - PrevH.m)) * ( dPda(t = age_vec, parameters = unname(fit_prev.m$coefficients)) + dPdt(a = age_vec, deriv = deriv_t)) +
      PrevH.m * predict(mort.spline.m, x = age_vec)$y
    inc.spline.m <- splinefun(x = age_vec, y = IncPrev.m)
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



                                             PrevH.f <- P(t = age_vec, parameters = unname(fit_prev.f$coefficients))
                                             IncPrev.f <- (1 / (1 - PrevH.f)) * ( dPda(t = age_vec, parameters = unname(fit_prev.f$coefficients))  + dPdt(a = age_vec, deriv = deriv_t)) +
                                               PrevH.f * predict(mort.spline.f, x = age_vec)$y

                                             inc.spline.f <- splinefun(x = age_vec, y = IncPrev.f)
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


                                             PrevH.m <- P(t = age_vec, parameters = unname(fit_prev.m$coefficients))
                                             IncPrev.m <- (1 / (1 - PrevH.m)) * ( dPda(t = age_vec, parameters = unname(fit_prev.m$coefficients)) + dPdt(a = age_vec, deriv = deriv_t)) +
                                               PrevH.m * predict(mort.spline.m, x = age_vec)$y

                                             inc.spline.m <- splinefun(x = age_vec, y = IncPrev.m)
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

averages.samp <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                            Age_L = c(15,20,25,30,15,25,15,15),
                            Age_U = c(20,25,30,35,25,35,30,35),
                            Women = rep(NA,8),
                            Women.LB = rep(NA,8),
                            Women.UB = rep(NA,8),
                            Men = rep(NA,8),
                            Men.LB = rep(NA,8),
                            Men.UB = rep(NA,8),
                            Both = rep(NA,8),
                            Both.LB = rep(NA,8),
                            Both.UB = rep(NA,8))
for (i in 1:nrow(averages.samp)) {
  Women <- average_incidence_sc(data = kzn,
                                sex = "Female",
                                age.range.fit = c(15,35),
                                step = 0.25,
                                age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                                weighting = "sampling",
                                poptable = poptable,
                                mort.spline.f = mort.spline.f,
                                mort.spline.m = mort.spline.m,
                                cores = cores,
                                n_bootstraps = n_bootstraps,
                                alpha = 0.05)
  Men <- average_incidence_sc(data = kzn,
                              sex = "Male",
                              age.range.fit = c(15,35),
                              step = 0.25,
                              age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                              weighting = "sampling",
                              poptable = poptable,
                              mort.spline.f = mort.spline.f,
                              mort.spline.m = mort.spline.m,
                              cores = cores,
                              n_bootstraps = n_bootstraps,
                              alpha = 0.05)
  Both <- average_incidence_sc(data = kzn,
                               sex = "both",
                               age.range.fit = c(15,35),
                               step = 0.25,
                               age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                               weighting = "sampling",
                               poptable = poptable,
                               mort.spline.f = mort.spline.f,
                               mort.spline.m = mort.spline.m,
                               cores = cores,
                               n_bootstraps = n_bootstraps,
                               alpha = 0.05)
  averages.samp$Women[i] <- Women[[1]]
  averages.samp$Women.LB[i] <- Women[[2]]
  averages.samp$Women.UB[i] <- Women[[3]]
  averages.samp$Men[i] <- Men[[1]]
  averages.samp$Men.LB[i] <- Men[[2]]
  averages.samp$Men.UB[i] <- Men[[3]]
  averages.samp$Both[i] <- Both[[1]]
  averages.samp$Both.LB[i] <- Both[[2]]
  averages.samp$Both.UB[i] <- Both[[3]]
}
write_feather(averages.samp, "averages.sc.sampweighted.feather")

# if(Sys.info()['login']=='eduardgrebe') {
#   average_incidence_sc(data = kzn,
#                     sex = "both",
#                     age.range.fit = c(15,35),
#                     step = 0.25,
#                     age.range.av = c(15,30),
#                     weighting = "population",
#                     poptable = poptable,
#                     mort.spline.f = mort.spline.f,
#                     mort.spline.m = mort.spline.m,
#                     cores = cores,
#                     n_bootstraps = n_bootstraps,
#                     alpha = 0.05,
#                     debug = TRUE)
# }

averages.pop <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                           Age_L = c(15,20,25,30,15,25,15,15),
                           Age_U = c(20,25,30,35,25,35,30,35),
                           Women = rep(NA,8),
                           Women.LB = rep(NA,8),
                           Women.UB = rep(NA,8),
                           Men = rep(NA,8),
                           Men.LB = rep(NA,8),
                           Men.UB = rep(NA,8),
                           Both = rep(NA,8),
                           Both.LB = rep(NA,8),
                           Both.UB = rep(NA,8))
for (i in 1:nrow(averages.pop)) {
  Women <- average_incidence_sc(data = kzn,
                                sex = "Female",
                                age.range.fit = c(15,35),
                                step = 0.25,
                                age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                                weighting = "population",
                                poptable = poptable,
                                mort.spline.f = mort.spline.f,
                                mort.spline.m = mort.spline.m,
                                cores = cores,
                                n_bootstraps = n_bootstraps,
                                alpha = 0.05)
  Men <- average_incidence_sc(data = kzn,
                              sex = "Male",
                              age.range.fit = c(15,35),
                              step = 0.25,
                              age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                              weighting = "population",
                              poptable = poptable,
                              mort.spline.f = mort.spline.f,
                              mort.spline.m = mort.spline.m,
                              cores = cores,
                              n_bootstraps = n_bootstraps,
                              alpha = 0.05)
  Both <- average_incidence_sc(data = kzn,
                               sex = "both",
                               age.range.fit = c(15,35),
                               step = 0.25,
                               age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                               weighting = "population",
                               poptable = poptable,
                               mort.spline.f = mort.spline.f,
                               mort.spline.m = mort.spline.m,
                               cores = cores,
                               n_bootstraps = n_bootstraps,
                               alpha = 0.05)
  averages.pop$Women[i] <- Women[[1]]
  averages.pop$Women.LB[i] <- Women[[2]]
  averages.pop$Women.UB[i] <- Women[[3]]
  averages.pop$Men[i] <- Men[[1]]
  averages.pop$Men.LB[i] <- Men[[2]]
  averages.pop$Men.UB[i] <- Men[[3]]
  averages.pop$Both[i] <- Both[[1]]
  averages.pop$Both.LB[i] <- Both[[2]]
  averages.pop$Both.UB[i] <- Both[[3]]
}
write_feather(averages.pop, "averages.sc.popweighted.feather")




## *** Age-specific biomarker method alone ***

average_incidence_bm <- function(data, sex = "Female", age.range.fit = c(15,35),
                                 step = 0.25, age.range.av = c(15,35), deriv_t = 0,
                                 weighting = "sampling", poptable = NULL, mdri, frr,
                                 bigT, cores = 4, n_bootstraps = 10000, alpha = 0.05,
                                 debug = FALSE)
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
    PrevR.f <- PR(t = age_vec, parameters = unname(fit_prevR.f$coefficients))
    IncR.f  <- IR(prevH = PrevH.f,
                  prevR = PrevR.f,
                  mdri = mdri,
                  frr = frr,
                  bigT = bigT)

    inc.spline.f <- splinefun(x = age_vec, y = IncR.f)
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
    PrevR.m <- PR(t = age_vec, parameters = unname(fit_prevR.m$coefficients))
    IncR.m  <- IR(prevH = PrevH.m,
                  prevR = PrevR.m,
                  mdri = mdri,
                  frr = frr,
                  bigT = bigT)

    inc.spline.m <- splinefun(x = age_vec, y = IncR.m)
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
                                             PrevR.f <- PR(t = age_vec, parameters = unname(fit_prevR.f$coefficients))
                                             IncR.f  <- IR(prevH = PrevH.f,
                                                           prevR = PrevR.f,
                                                           mdri = mdri,
                                                           frr = frr,
                                                           bigT = bigT)

                                             inc.spline.f <- splinefun(x = age_vec, y = IncR.f)
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
                                             PrevR.m <- PR(t = age_vec, parameters = unname(fit_prevR.m$coefficients))
                                             IncR.m  <- IR(prevH = PrevH.m,
                                                           prevR = PrevR.m,
                                                           mdri = mdri,
                                                           frr = frr,
                                                           bigT = bigT)

                                             inc.spline.m <- splinefun(x = age_vec, y = IncR.m)
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

averages.samp <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                            Age_L = c(15,20,25,30,15,25,15,15),
                            Age_U = c(20,25,30,35,25,35,30,35),
                            Women = rep(NA,8),
                            Women.LB = rep(NA,8),
                            Women.UB = rep(NA,8),
                            Men = rep(NA,8),
                            Men.LB = rep(NA,8),
                            Men.UB = rep(NA,8),
                            Both = rep(NA,8),
                            Both.LB = rep(NA,8),
                            Both.UB = rep(NA,8))
for (i in 1:nrow(averages.samp)) {
  Women <- average_incidence_bm(data = kzn,
                                sex = "Female",
                                age.range.fit = c(15,35),
                                step = 0.25,
                                age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                                weighting = "sampling",
                                poptable = poptable,
                                mdri = mdriYr,
                                frr = frrContext,
                                bigT = 2,
                                cores = cores,
                                n_bootstraps = n_bootstraps,
                                alpha = 0.05)
  Men <- average_incidence_bm(data = kzn,
                              sex = "Male",
                              age.range.fit = c(15,35),
                              step = 0.25,
                              age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                              weighting = "sampling",
                              poptable = poptable,
                              mdri = mdriYr,
                              frr = frrContext,
                              bigT = 2,
                              cores = cores,
                              n_bootstraps = n_bootstraps,
                              alpha = 0.05)
  Both <- average_incidence_bm(data = kzn,
                               sex = "both",
                               age.range.fit = c(15,35),
                               step = 0.25,
                               age.range.av = c(averages.samp$Age_L[i],averages.samp$Age_U[i]),
                               weighting = "sampling",
                               poptable = poptable,
                               mdri = mdriYr,
                               frr = frrContext,
                               bigT = 2,
                               cores = cores,
                               n_bootstraps = n_bootstraps,
                               alpha = 0.05)
  averages.samp$Women[i] <- Women[[1]]
  averages.samp$Women.LB[i] <- Women[[2]]
  averages.samp$Women.UB[i] <- Women[[3]]
  averages.samp$Men[i] <- Men[[1]]
  averages.samp$Men.LB[i] <- Men[[2]]
  averages.samp$Men.UB[i] <- Men[[3]]
  averages.samp$Both[i] <- Both[[1]]
  averages.samp$Both.LB[i] <- Both[[2]]
  averages.samp$Both.UB[i] <- Both[[3]]
}
write_feather(averages.samp, "averages.bm.sampweighted.feather")

averages.pop <- data_frame(AgeRange = c("[15,20)","[20,25)","[25,30)","[30,35)","[15,25)","[25,35)","[15,30)","[15,35)"),
                           Age_L = c(15,20,25,30,15,25,15,15),
                           Age_U = c(20,25,30,35,25,35,30,35),
                           Women = rep(NA,8),
                           Women.LB = rep(NA,8),
                           Women.UB = rep(NA,8),
                           Men = rep(NA,8),
                           Men.LB = rep(NA,8),
                           Men.UB = rep(NA,8),
                           Both = rep(NA,8),
                           Both.LB = rep(NA,8),
                           Both.UB = rep(NA,8))
for (i in 1:nrow(averages.pop)) {
  Women <- average_incidence_bm(data = kzn,
                                sex = "Female",
                                age.range.fit = c(15,35),
                                step = 0.25,
                                age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                                weighting = "population",
                                poptable = poptable,
                                mdri = mdriYr,
                                frr = frrContext,
                                bigT = 2,
                                cores = cores,
                                n_bootstraps = n_bootstraps,
                                alpha = 0.05)
  Men <- average_incidence_bm(data = kzn,
                              sex = "Male",
                              age.range.fit = c(15,35),
                              step = 0.25,
                              age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                              weighting = "population",
                              poptable = poptable,
                              mdri = mdriYr,
                              frr = frrContext,
                              bigT = 2,
                              cores = cores,
                              n_bootstraps = n_bootstraps,
                              alpha = 0.05)
  Both <- average_incidence_bm(data = kzn,
                               sex = "both",
                               age.range.fit = c(15,35),
                               step = 0.25,
                               age.range.av = c(averages.pop$Age_L[i],averages.pop$Age_U[i]),
                               weighting = "population",
                               poptable = poptable,
                               mdri = mdriYr,
                               frr = frrContext,
                               bigT = 2,
                               cores = cores,
                               n_bootstraps = n_bootstraps,
                               alpha = 0.05)
  averages.pop$Women[i] <- Women[[1]]
  averages.pop$Women.LB[i] <- Women[[2]]
  averages.pop$Women.UB[i] <- Women[[3]]
  averages.pop$Men[i] <- Men[[1]]
  averages.pop$Men.LB[i] <- Men[[2]]
  averages.pop$Men.UB[i] <- Men[[3]]
  averages.pop$Both[i] <- Both[[1]]
  averages.pop$Both.LB[i] <- Both[[2]]
  averages.pop$Both.UB[i] <- Both[[3]]
}
write_feather(averages.pop, "averages.bm.popweighted.feather")
