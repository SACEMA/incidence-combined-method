# Copyright (C) 2018 Eduard Grebe, Alex Welte and Stellenbosch University
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

 setwd("~/dev/incidence-kzn/")
  library(feather); library(haven); library(dplyr); library(inctools);
  library(glm2); library(parallel); library(doParallel); library(foreach);
  library(readxl); library(stringr); library(lubridate); library(tidyr);
  library(survey); library(import); library(abind);

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

mort.spline.zero <- smooth.spline(x = seq(15,60,1), y = rep(0,46))

poptable <- read_feather("poptable.feather")

load("bs_params_w1535.RData")
load("bs_params_m1535.RData")
load("bs_params_b1535.RData")


### *** FRR sensitivity analysis *** ###
  n_bootstraps <- 10000 # 100
  cores <- 4

# Standard parameters

mdri <- 224.5231 - 7.71 # Adjustment for NAT screening with DT 100 c/ml, pools of 5 (effective DT 500c/ml) # Subtype C!!
mdriYr <- mdri / 365.25
rse_mdri <- 0.0619132331026123 # NB: uncertainty in NAT adjustment not known or taken into account
frrContext <- 0.00166780137766038
rse_frrContext <- 0.465416630819573

FRRs <- seq(0,0.01,by = 0.001)
RSE_FRR <- 0.5

incidences_frr_sens <- read_feather("incidences_frr_sens_empty.feather")
for (frr_iterate in FRRs) {
  # Women 15-35
 incarray_w1535 <- extract_incidences(params = bs_params_w1535, 
                                      lambda = 0, 
                                      mort.spline = mort.spline.f, 
                                      age.range = c(15,35), 
                                      age.step = 1, 
                                      cores = cores, 
                                      mdri = mdriYr, 
                                      frr = frr_iterate, 
                                      bigT = 2)
 inc_w1535 <- incidence(data = kzn, 
                           sex = "Female", 
                           age.range = c(15,35), 
                           incmat = incarray_w1535, 
                           lambda = 0, 
                           mort.spline = mort.spline.f,
                           mdri = mdriYr, 
                           rse_mdri = rse_mdri,
                           frr = frr_iterate, 
                           rse_frr = RSE_FRR,
                           bigT = 2, 
                           alpha = 0.05, 
                        debug = FALSE)
  
  inc_w1535$IncPrev.LB <- inc_w1535$IncPrev - qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  inc_w1535$IncPrev.LB <- ifelse(inc_w1535$IncPrev.LB<0,0,inc_w1535$IncPrev.LB)
  inc_w1535$IncPrev.UB <- inc_w1535$IncPrev + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$IncR.LB <- inc_w1535$IncR - qnorm(1-(0.05/2))*inc_w1535$IncR.sigma
  inc_w1535$IncR.LB <- ifelse(inc_w1535$IncR.LB<0,0,inc_w1535$IncR.LB)
  inc_w1535$IncR.UB <- inc_w1535$IncR + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$Sex <- "Women"
  inc_w1535$FRR <- frr_iterate
  
  #save(incarray_w1535, file = paste0("incarray_w1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_w1535, paste0("inc_w1535_FRR_",frr_iterate,".feather"))
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_w1535)  
  
  # Men 15-35
  incarray_m1535 <- extract_incidences(params = bs_params_m1535, 
                                       lambda = 0, 
                                       mort.spline = mort.spline.m, 
                                       age.range = c(15,35), 
                                       age.step = 1, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frr_iterate, 
                                       bigT = 2)
  inc_m1535 <- incidence(data = kzn, 
                         sex = "Male", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_m1535, 
                         lambda = 0, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frr_iterate, 
                         rse_frr = RSE_FRR,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_m1535$IncPrev.LB <- inc_m1535$IncPrev - qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  inc_m1535$IncPrev.LB <- ifelse(inc_m1535$IncPrev.LB<0,0,inc_m1535$IncPrev.LB)
  inc_m1535$IncPrev.UB <- inc_m1535$IncPrev + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$IncR.LB <- inc_m1535$IncR - qnorm(1-(0.05/2))*inc_m1535$IncR.sigma
  inc_m1535$IncR.LB <- ifelse(inc_m1535$IncR.LB<0,0,inc_m1535$IncR.LB)
  inc_m1535$IncR.UB <- inc_m1535$IncR + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$Sex <- "Men"
  inc_m1535$FRR <- frr_iterate
  
  #save(incarray_m1535, file = paste0("incarray_m1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_m1535, paste0("inc_m1535_FRR_",frr_iterate,".feather"))
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_m1535)
  
  
  # All 15-35
  incarray_b1535 <- extract_incidences(params = bs_params_b1535, 
                                       lambda = 0, 
                                       mort.spline = mort.spline.b, 
                                       age.range = c(15,35), 
                                       age.step = 1, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frr_iterate, 
                                       bigT = 2)
  inc_b1535 <- incidence(data = kzn, 
                         sex = "both", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_b1535, 
                         lambda = 0, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frr_iterate, 
                         rse_frr = RSE_FRR,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_b1535$IncPrev.LB <- inc_b1535$IncPrev - qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  inc_b1535$IncPrev.LB <- ifelse(inc_b1535$IncPrev.LB<0,0,inc_b1535$IncPrev.LB)
  inc_b1535$IncPrev.UB <- inc_b1535$IncPrev + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$IncR.LB <- inc_b1535$IncR - qnorm(1-(0.05/2))*inc_b1535$IncR.sigma
  inc_b1535$IncR.LB <- ifelse(inc_b1535$IncR.LB<0,0,inc_b1535$IncR.LB)
  inc_b1535$IncR.UB <- inc_b1535$IncR + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$Sex <- "Both"
  inc_b1535$FRR <- frr_iterate
 
  #save(incarray_b1535, file = paste0("incarray_b1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_b1535, paste0("inc_b1535_FRR_",frr_iterate,".feather"))
  
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_b1535)
}

write_feather(incidences_frr_sens, "incidences_frr_sens_rse1.feather")


FRRs <- seq(0,0.01,by = 0.001)
RSE_FRR <- 0

incidences_frr_sens <- read_feather("incidences_frr_sens_empty.feather")
for (frr_iterate in FRRs) {
  # Women 15-35
  incarray_w1535 <- extract_incidences(params = bs_params_w1535, 
                                       lambda = 0, 
                                       mort.spline = mort.spline.f, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frr_iterate, 
                                       bigT = 2)
  inc_w1535 <- incidence(data = kzn, 
                         sex = "Female", 
                         age.range = c(15,35), 
                         incmat = incarray_w1535, 
                         lambda = 0, 
                         mort.spline = mort.spline.f,
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frr_iterate, 
                         rse_frr = RSE_FRR,
                         bigT = 2, 
                         alpha = 0.05, 
                         debug = FALSE)
  
  inc_w1535$IncPrev.LB <- inc_w1535$IncPrev - qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  inc_w1535$IncPrev.LB <- ifelse(inc_w1535$IncPrev.LB<0,0,inc_w1535$IncPrev.LB)
  inc_w1535$IncPrev.UB <- inc_w1535$IncPrev + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$IncR.LB <- inc_w1535$IncR - qnorm(1-(0.05/2))*inc_w1535$IncR.sigma
  inc_w1535$IncR.LB <- ifelse(inc_w1535$IncR.LB<0,0,inc_w1535$IncR.LB)
  inc_w1535$IncR.UB <- inc_w1535$IncR + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$Sex <- "Women"
  inc_w1535$FRR <- frr_iterate
  
  #save(incarray_w1535, file = paste0("incarray_w1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_w1535, paste0("inc_w1535_FRR_",frr_iterate,".feather"))
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_w1535)  
  
  # Men 15-35
  incarray_m1535 <- extract_incidences(params = bs_params_m1535, 
                                       lambda = 0, 
                                       mort.spline = mort.spline.m, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frr_iterate, 
                                       bigT = 2)
  inc_m1535 <- incidence(data = kzn, 
                         sex = "Male", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_m1535, 
                         lambda = 0, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frr_iterate, 
                         rse_frr = RSE_FRR,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_m1535$IncPrev.LB <- inc_m1535$IncPrev - qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  inc_m1535$IncPrev.LB <- ifelse(inc_m1535$IncPrev.LB<0,0,inc_m1535$IncPrev.LB)
  inc_m1535$IncPrev.UB <- inc_m1535$IncPrev + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$IncR.LB <- inc_m1535$IncR - qnorm(1-(0.05/2))*inc_m1535$IncR.sigma
  inc_m1535$IncR.LB <- ifelse(inc_m1535$IncR.LB<0,0,inc_m1535$IncR.LB)
  inc_m1535$IncR.UB <- inc_m1535$IncR + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$Sex <- "Men"
  inc_m1535$FRR <- frr_iterate
  
  #save(incarray_m1535, file = paste0("incarray_m1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_m1535, paste0("inc_m1535_FRR_",frr_iterate,".feather"))
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_m1535)
  
  
  # All 15-35
  incarray_b1535 <- extract_incidences(params = bs_params_b1535, 
                                       lambda = 0, 
                                       mort.spline = mort.spline.b, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frr_iterate, 
                                       bigT = 2)
  inc_b1535 <- incidence(data = kzn, 
                         sex = "both", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_b1535, 
                         lambda = 0, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frr_iterate, 
                         rse_frr = RSE_FRR,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_b1535$IncPrev.LB <- inc_b1535$IncPrev - qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  inc_b1535$IncPrev.LB <- ifelse(inc_b1535$IncPrev.LB<0,0,inc_b1535$IncPrev.LB)
  inc_b1535$IncPrev.UB <- inc_b1535$IncPrev + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$IncR.LB <- inc_b1535$IncR - qnorm(1-(0.05/2))*inc_b1535$IncR.sigma
  inc_b1535$IncR.LB <- ifelse(inc_b1535$IncR.LB<0,0,inc_b1535$IncR.LB)
  inc_b1535$IncR.UB <- inc_b1535$IncR + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$Sex <- "Both"
  inc_b1535$FRR <- frr_iterate
  
  #save(incarray_b1535, file = paste0("incarray_b1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_b1535, paste0("inc_b1535_FRR_",frr_iterate,".feather"))
  
  incidences_frr_sens <- bind_rows(incidences_frr_sens,inc_b1535)
}

write_feather(incidences_frr_sens, "incidences_frr_sens_rse0.feather")




### *** dPdt sensitivity analysis *** ###

halflife <- function(lambda) {
  return(-log(0.5)/lambda)
}

incidences_dpdt_sens <- read_feather("incidences_frr_sens_empty.feather")
incidences_dpdt_sens <- incidences_dpdt_sens %>%
  rename(lambda = FRR)

# Standard params

mdri <- 224.5231 - 7.71 # Adjustment for NAT screening with DT 100 c/ml, pools of 5 (effective DT 500c/ml) # Subtype C!!
mdriYr <- mdri / 365.25
rse_mdri <- 0.0619132331026123 # NB: uncertainty in NAT adjustment not known or taken into account
frrContext <- 0.00166780137766038
rse_frrContext <- 0.465416630819573

lambdas <- c(-0.07,-0.035,0.035,0.07,0.14)

for (lambda_iterate in lambdas) {
  # Women 15-35
  incarray_w1535 <- extract_incidences(params = bs_params_w1535, 
                                       lambda = lambda_iterate, 
                                       mort.spline = mort.spline.f, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frrContext, 
                                       bigT = 2)
  inc_w1535 <- incidence(data = kzn, 
                         sex = "Female", 
                         age.range = c(15,35), 
                         incmat = incarray_w1535, 
                         lambda = lambda_iterate, 
                         mort.spline = mort.spline.f,
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frrContext, 
                         rse_frr = rse_frrContext,
                         bigT = 2, 
                         alpha = 0.05, 
                         debug = FALSE)
  
  inc_w1535$IncPrev.LB <- inc_w1535$IncPrev - qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  inc_w1535$IncPrev.LB <- ifelse(inc_w1535$IncPrev.LB<0,0,inc_w1535$IncPrev.LB)
  inc_w1535$IncPrev.UB <- inc_w1535$IncPrev + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$IncR.LB <- inc_w1535$IncR - qnorm(1-(0.05/2))*inc_w1535$IncR.sigma
  inc_w1535$IncR.LB <- ifelse(inc_w1535$IncR.LB<0,0,inc_w1535$IncR.LB)
  inc_w1535$IncR.UB <- inc_w1535$IncR + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
  
  inc_w1535$Sex <- "Women"
  inc_w1535$lambda <- lambda_iterate
  
  #save(incarray_w1535, file = paste0("incarray_w1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_w1535, paste0("inc_w1535_FRR_",frr_iterate,".feather"))
  incidences_dpdt_sens <- bind_rows(incidences_dpdt_sens,inc_w1535)  
  
  # Men 15-35
  incarray_m1535 <- extract_incidences(params = bs_params_m1535, 
                                       lambda = lambda_iterate, 
                                       mort.spline = mort.spline.m, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frrContext, 
                                       bigT = 2)
  inc_m1535 <- incidence(data = kzn, 
                         sex = "Male", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_m1535, 
                         lambda = lambda_iterate, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frrContext, 
                         rse_frr = rse_frrContext,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_m1535$IncPrev.LB <- inc_m1535$IncPrev - qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  inc_m1535$IncPrev.LB <- ifelse(inc_m1535$IncPrev.LB<0,0,inc_m1535$IncPrev.LB)
  inc_m1535$IncPrev.UB <- inc_m1535$IncPrev + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$IncR.LB <- inc_m1535$IncR - qnorm(1-(0.05/2))*inc_m1535$IncR.sigma
  inc_m1535$IncR.LB <- ifelse(inc_m1535$IncR.LB<0,0,inc_m1535$IncR.LB)
  inc_m1535$IncR.UB <- inc_m1535$IncR + qnorm(1-(0.05/2))*inc_m1535$IncPrev.sigma
  
  inc_m1535$Sex <- "Men"
  inc_m1535$lambda <- lambda_iterate
  
  #save(incarray_m1535, file = paste0("incarray_m1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_m1535, paste0("inc_m1535_FRR_",frr_iterate,".feather"))
  incidences_dpdt_sens <- bind_rows(incidences_dpdt_sens,inc_m1535)
  
  
  # All 15-35
  incarray_b1535 <- extract_incidences(params = bs_params_b1535, 
                                       lambda = lambda_iterate, 
                                       mort.spline = mort.spline.b, 
                                       age.range = c(15,35), 
                                       age.step = 5, 
                                       cores = cores, 
                                       mdri = mdriYr, 
                                       frr = frrContext, 
                                       bigT = 2)
  inc_b1535 <- incidence(data = kzn, 
                         sex = "both", 
                         age.range = c(15,35), 
                         mort.spline = mort.spline.m,
                         incmat = incarray_b1535, 
                         lambda = lambda_iterate, 
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frrContext, 
                         rse_frr = rse_frrContext,
                         bigT = 2, 
                         alpha = 0.05)
  
  inc_b1535$IncPrev.LB <- inc_b1535$IncPrev - qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  inc_b1535$IncPrev.LB <- ifelse(inc_b1535$IncPrev.LB<0,0,inc_b1535$IncPrev.LB)
  inc_b1535$IncPrev.UB <- inc_b1535$IncPrev + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$IncR.LB <- inc_b1535$IncR - qnorm(1-(0.05/2))*inc_b1535$IncR.sigma
  inc_b1535$IncR.LB <- ifelse(inc_b1535$IncR.LB<0,0,inc_b1535$IncR.LB)
  inc_b1535$IncR.UB <- inc_b1535$IncR + qnorm(1-(0.05/2))*inc_b1535$IncPrev.sigma
  
  inc_b1535$Sex <- "Both"
  inc_b1535$lambda <- lambda_iterate
  
  #save(incarray_b1535, file = paste0("incarray_b1535_FRR_",frr_iterate,".RData"))
  #write_feather(inc_b1535, paste0("inc_b1535_FRR_",frr_iterate,".feather"))
  
  incidences_dpdt_sens <- bind_rows(incidences_dpdt_sens,inc_b1535)
}

write_feather(incidences_dpdt_sens, "incidences_dpdt_sens.feather")


