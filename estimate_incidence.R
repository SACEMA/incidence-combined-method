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

  cores <- 4
  
kzn <- read_feather("kzn.feather")


mortality_table <- read_feather("mortality_table.feather")

mortality_f_2013 <- mortality_table %>%
  filter(Sex == "Female") %>%
  select(age = Age, mort = `2013`)
mort.spline.f <- smooth.spline(x = mortality_f_2013$age, y = mortality_f_2013$mort)

mortality_m_2013 <- mortality_table %>%
  filter(Sex == "Male") %>%
  select(age = Age, mort = `2013`)
mort.spline.m <- smooth.spline(x = mortality_m_2013$age, y = mortality_m_2013$mort)

mortality_b_2013 <- mortality_table %>%
  filter(Sex == "Combined") %>%
  select(age = Age, mort = `2013`)
mort.spline.b <- smooth.spline(x = mortality_b_2013$age, y = mortality_b_2013$mort)

mort.spline.zero <- smooth.spline(x = seq(15,60,1), y = rep(0,46))


# Parameters

# Use findings from earlier calibration:
# recent: NAT yield | (LAg < 1.75 & BioRad < 20 & VL > 75)
mdri <- 224.5231 - 7.71 # Adjustment for NAT screening with DT 100 c/ml, pools of 5 (effective DT 500c/ml) # Subtype C!!
mdriYr <- mdri / 365.25
rse_mdri <- 0.0619132331026123 # NB: uncertainty in NAT adjustment not known or taken into account
frrContext <- 0.00166780137766038
rse_frrContext <- 0.465416630819573

# Women 15-35
  bs_params_w1535 <- bootstrap_params(data = kzn, sex = "Female", age.range = c(15,35), n_bootstraps = 12500, cores = cores)
  incarray_w1535 <- extract_incidences(params = bs_params_w1535, deriv_t = 0, 
                                       mort.spline = mort.spline.f, 
                                       age.range = c(15,35), age.step = 0.25, 
                                       cores = cores, mdri = mdriYr, 
                                       frr = frrContext, bigT = 2)

  inc_w1535 <- incidence(data = kzn, 
                         sex = "Female", 
                         age.range = c(15,35), 
                         incmat = incarray_w1535, 
                         deriv_t = 0, 
                         mort.spline = mort.spline.f,
                         mdri = mdriYr, 
                         rse_mdri = rse_mdri,
                         frr = frrContext, 
                         rse_frr = rse_frrContext,
                         bigT = 2, 
                         alpha = 0.05, debug = FALSE)




inc_w1535$IncPrev.LB <- inc_w1535$IncPrev - qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma
inc_w1535$IncPrev.LB <- ifelse(inc_w1535$IncPrev.LB<0,0,inc_w1535$IncPrev.LB)
inc_w1535$IncPrev.UB <- inc_w1535$IncPrev + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma

inc_w1535$IncR.LB <- inc_w1535$IncR - qnorm(1-(0.05/2))*inc_w1535$IncR.sigma
inc_w1535$IncR.LB <- ifelse(inc_w1535$IncR.LB<0,0,inc_w1535$IncR.LB)
inc_w1535$IncR.UB <- inc_w1535$IncR + qnorm(1-(0.05/2))*inc_w1535$IncPrev.sigma

inc_w1535_long <- data_frame(Age = rep(inc_w1535$Age,3), 
                             Estimate = factor(c(rep("Incidence (synthetic cohort)",nrow(inc_w1535)),
                                                 rep("Incidence (biomarker)",nrow(inc_w1535)),
                                                 rep("Weighted average",nrow(inc_w1535)))), 
                             Incidence = c(inc_w1535$IncPrev,
                                           inc_w1535$IncR,
                                           inc_w1535$IncW), 
                             LB = c(inc_w1535$IncPrev.LB,
                                    inc_w1535$IncR.LB,
                                    inc_w1535$IncW.LB), 
                             UB = c(inc_w1535$IncPrev.UB,
                                    inc_w1535$IncR.UB,
                                    inc_w1535$IncW.UB)
)

women_summary <- inc_w1535 %>%
  filter(Age %in% seq(15,35,1))


save(bs_params_w1535, file = "bs_params_w1535.RData")
save(incarray_w1535, file = "incarray_w1535.RData")
write_feather(inc_w1535, "inc_w1535.feather")
write_feather(women_summary, "women_summary.feather")



# Men 15-35
bs_params_m1535 <- bootstrap_params(data = kzn, sex = "Male", age.range = c(15,35), n_bootstraps = 12500, cores = cores)
incarray_m1535 <- extract_incidences(params = bs_params_m1535, deriv_t = 0, mort.spline = mort.spline.m, age.range = c(15,35), age.step = 0.25, cores = cores, mdri = mdriYr, frr = frrContext, bigT = 2)
inc_m1535 <- incidence(data = kzn, 
                       sex = "Male", 
                       age.range = c(15,35), 
                       mort.spline = mort.spline.m,
                       incmat = incarray_m1535, 
                       deriv_t = 0, 
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

inc_m1535_long <- data_frame(Age = rep(inc_m1535$Age,3), 
                             Estimate = factor(c(rep("Incidence (synthetic cohort)",nrow(inc_m1535)),
                                                 rep("Incidence (biomarker)",nrow(inc_m1535)),
                                                 rep("Weighted average",nrow(inc_m1535)))), 
                             Incidence = c(inc_m1535$IncPrev,
                                           inc_m1535$IncR,
                                           inc_m1535$IncW), 
                             LB = c(inc_m1535$IncPrev.LB,
                                    inc_m1535$IncR.LB,
                                    inc_m1535$IncW.LB), 
                             UB = c(inc_m1535$IncPrev.UB,
                                    inc_m1535$IncR.UB,
                                    inc_m1535$IncW.UB)
)

men_summary <- inc_m1535 %>%
  filter(Age %in% seq(15,35,1))

save(bs_params_m1535, file = "bs_params_m1535.RData")
save(incarray_m1535, file = "incarray_m1535.RData")
write_feather(inc_m1535, "inc_m1535.feather")
write_feather(men_summary, "men_summary.feather")


# All 15-35
bs_params_b1535 <- bootstrap_params(data = kzn, sex = "both", age.range = c(15,35), n_bootstraps = 12500, cores = cores)
incarray_b1535 <- extract_incidences(params = bs_params_b1535, deriv_t = 0, mort.spline = mort.spline.b, age.range = c(15,35), age.step = 0.25, cores = cores, mdri = mdriYr, frr = frrContext, bigT = 2)
inc_b1535 <- incidence(data = kzn, 
                       sex = "both", 
                       age.range = c(15,35), 
                       mort.spline = mort.spline.m,
                       incmat = incarray_b1535, 
                       deriv_t = 0, 
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

inc_b1535_long <- data_frame(Age = rep(inc_b1535$Age,3), 
                             Estimate = factor(c(rep("Incidence (synthetic cohort)",nrow(inc_b1535)),
                                                 rep("Incidence (biomarker)",nrow(inc_b1535)),
                                                 rep("Weighted average",nrow(inc_b1535)))), 
                             Incidence = c(inc_b1535$IncPrev,
                                           inc_b1535$IncR,
                                           inc_b1535$IncW), 
                             LB = c(inc_b1535$IncPrev.LB,
                                    inc_b1535$IncR.LB,
                                    inc_b1535$IncW.LB), 
                             UB = c(inc_b1535$IncPrev.UB,
                                    inc_b1535$IncR.UB,
                                    inc_b1535$IncW.UB)
)

both_summary <- inc_b1535 %>%
  filter(Age %in% seq(15,35,1))

save(bs_params_b1535, file = "bs_params_b1535.RData")
save(incarray_b1535, file = "incarray_b1535.RData")
write_feather(inc_b1535, "inc_b1535.feather")
write_feather(both_summary, "both_summary.feather")