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

knitr::opts_knit$set(root.dir = "~/dev/incidence-kzn/")
library(haven)
library(dplyr)
library(inctools)
library(glm2)
library(parallel)
library(doParallel)
library(foreach)
library(readxl)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(survey)
library(feather)
library(readr)


###*** Primary results ***###

load("output/bs_params_w1535.RData")
load("output/incarray_w1535.RData")
inc_w1535 <- read_feather("output/inc_w1535.feather")
women_summary <- read_feather("output/women_summary.feather")
women_summary$Sex <- "Women"

load("output/bs_params_m1535.RData")
load("output/incarray_m1535.RData")
inc_m1535 <- read_feather("output/inc_m1535.feather")
men_summary <- read_feather("output/men_summary.feather")
men_summary$Sex <- "Men"

load("output/bs_params_b1535.RData")
load("output/incarray_b1535.RData")
inc_b1535 <- read_feather("output/inc_b1535.feather")
both_summary <- read_feather("output/both_summary.feather")
both_summary$Sex <- "Both"

summary <- bind_rows(men_summary, women_summary, both_summary)
write_delim(summary, "output/incidence_summary.csv", delim = ";")

inc_w1535_long <- data_frame(Age = rep(inc_w1535$Age,3),
                             Method = factor(c(rep("Synthetic cohort",nrow(inc_w1535)),
                                                 rep("Biomarker",nrow(inc_w1535)),
                                                 rep("Combined",nrow(inc_w1535))),
                                               levels = c("Biomarker","Synthetic cohort","Combined")),
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


plt_w <- ggplot() +
  geom_ribbon(data = inc_w1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,30), breaks = seq(0,30,2)) +
  scale_x_continuous(limits = c(15,35), breaks = seq(15,35,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_w)
ggsave("output/Incidence_women_15-35_compare.pdf", w=7,h=4.5, dpi=2400)


plt_w <- ggplot() +
  geom_ribbon(data = inc_w1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_w1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) +
  scale_x_continuous(limits = c(15,30), breaks = seq(15,30,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-30\nFit to data for ages 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_w)
ggsave("output/Figure3.pdf", w=7,h=4.5, dpi=2400)
ggsave("output/Figure3.TIF", device = "tiff", w=14, h=9, units = "cm", dpi=600)

inc_m1535_long <- data_frame(Age = rep(inc_m1535$Age,3),
                             Method = factor(c(rep("Synthetic cohort",nrow(inc_m1535)),
                                                 rep("Biomarker",nrow(inc_m1535)),
                                                 rep("Combined",nrow(inc_m1535))),
                                               levels = c("Biomarker","Synthetic cohort","Combined")),
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

plt_m <- ggplot() +
  geom_ribbon(data = inc_m1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,15), breaks = seq(0,15,1)) +
  scale_x_continuous(limits = c(15,35), breaks = seq(15,35,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_m)
ggsave("output/Incidence_men_15-35_compare.pdf", w=7,h=4.5, dpi=2400)


plt_m <- ggplot() +
  geom_ribbon(data = inc_m1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_m1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) +
  scale_x_continuous(limits = c(15,30), breaks = seq(15,30,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-30\nFit to data for ages 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_m)
ggsave("output/Figure2.pdf", w=7,h=4.5, dpi=2400)
ggsave("output/Figure2.TIF", device = "tiff", w=14, h=9, units = "cm", dpi=600)


inc_b1535_long <- data_frame(Age = rep(inc_b1535$Age,3),
                             Method = factor(c(rep("Synthetic cohort",nrow(inc_b1535)),
                                                 rep("Biomarker",nrow(inc_b1535)),
                                                 rep("Combined",nrow(inc_b1535))),
                                               levels = c("Biomarker","Synthetic cohort","Combined")),
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

plt_b <- ggplot() +
  geom_ribbon(data = inc_b1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,18.5), breaks = seq(0,18,2)) +
  scale_x_continuous(limits = c(15,35), breaks = seq(15,35,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_b)
ggsave("output/Figure4.pdf", w=7,h=4.5, dpi=2400)
ggsave("output/Figure4.TIF", device = "tiff", w=14, h=9, units = "cm", dpi=600)

plt_b <- ggplot() +
  geom_ribbon(data = inc_b1535_long, aes(x = Age, ymin = LB*100, ymax = UB*100, fill = Method), alpha = 0.2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y = LB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y = UB*100, colour = Method), linetype = 2) +
  geom_line(data = inc_b1535_long, aes(x = Age, y=Incidence*100, colour = Method), size=0.8) +

  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) +
  scale_x_continuous(limits = c(15,30), breaks = seq(15,30,1)) +
  xlab("Age") +
  ylab("Incidence (cases/100PY)") +
  theme_bw() +
  #ggtitle("Incidence in women aged 15-30\nFit to data for ages 15-35") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "bottom")
print(plt_b)
ggsave("output/Figure1.pdf", w=7,h=4.5, dpi=2400)
ggsave("output/Figure1.TIF", device = "tiff", w=14, h=9, units = "cm", dpi=600)






###*** Sensitivity analyses - FRR ***###

# Biomarker-based estimates

import::from(gridExtra, grid.arrange)
sens_frr <- read_feather("output/incidences_frr_sens_rse1.feather")
write_delim(sens_frr, "output/incidences_frr_sens_rse1.csv", delim = ";")

sens_frr_plot <- sens_frr %>%
  select(Age, Sex, FRR, IncR, IncR.LB, IncR.UB) %>%
  mutate(FRR = FRR*100,
         IncR = IncR*100,
         IncR.LB = IncR.LB*100,
         IncR.UB = IncR.UB*100) %>%
  filter(Age %in% c(15,20,25,30))


sens_frr_plot <- sens_frr_plot %>%
  filter(Sex != "Both") %>%
  mutate(Sex = factor(case_when(
    Sex == "Men" ~ "Males",
    Sex == "Women" ~ "Females"
  ), levels = c("Males","Females")))

plot_15 <- ggplot(data = filter(sens_frr_plot, Age == 15), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncR.LB, ymax = IncR.UB, colour = Age
  #geom_errorbar(aes(ymin = IncR.LB, ymax = IncR.UB)) +
  geom_point(aes(y = IncR)) +
  geom_smooth(aes(y = IncR), se = F) +
  geom_smooth(aes(y = IncR.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncR.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("A: Males and females aged 15") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_15)

plot_20 <- ggplot(data = filter(sens_frr_plot, Age == 20), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncR.LB, ymax = IncR.UB, colour = Age
  #geom_errorbar(aes(ymin = IncR.LB, ymax = IncR.UB)) +
  geom_point(aes(y = IncR)) +
  geom_smooth(aes(y = IncR), se = F) +
  geom_smooth(aes(y = IncR.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncR.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("B: Males and females aged 20") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_20)

plot_25 <- ggplot(data = filter(sens_frr_plot, Age == 25), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncR.LB, ymax = IncR.UB, colour = Age
  #geom_errorbar(aes(ymin = IncR.LB, ymax = IncR.UB)) +
  geom_point(aes(y = IncR)) +
  geom_smooth(aes(y = IncR), se = F) +
  geom_smooth(aes(y = IncR.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncR.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("C: Males and females aged 25") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_25)

plot_30 <- ggplot(data = filter(sens_frr_plot, Age == 30), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncR.LB, ymax = IncR.UB, colour = Age
  #geom_errorbar(aes(ymin = IncR.LB, ymax = IncR.UB)) +
  geom_point(aes(y = IncR)) +
  geom_smooth(aes(y = IncR), se = F) +
  geom_smooth(aes(y = IncR.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncR.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("D: Males and females aged 30") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_30)

plot_biomarker <- grid.arrange(plot_15, plot_20, plot_25, plot_30, ncol = 2, nrow = 2)
ggsave("output/FigureS2_1.pdf", plot = plot_biomarker, w=14,h=16, dpi=2400)


# Combined method

# RSE on FRR 50%

import::from(gridExtra, grid.arrange)
sens_frr <- read_feather("output/incidences_frr_sens_rse1.feather")

sens_frr_plot <- sens_frr %>%
  select(Age, Sex, FRR, IncW, IncW.LB, IncW.UB) %>%
  mutate(FRR = FRR*100,
         IncW = IncW*100,
         IncW.LB = IncW.LB*100,
         IncW.UB = IncW.UB*100) %>%
  filter(Age %in% c(15,20,25,30))


sens_frr_plot <- sens_frr_plot %>%
  filter(Sex != "Both") %>%
  mutate(Sex = factor(case_when(
    Sex == "Men" ~ "Males",
    Sex == "Women" ~ "Females"
  ), levels = c("Males","Females")))

plot_15 <- ggplot(data = filter(sens_frr_plot, Age == 15), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("A: Males and females aged 15") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_15)

plot_20 <- ggplot(data = filter(sens_frr_plot, Age == 20), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("B: Males and females aged 20") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_20)

plot_25 <- ggplot(data = filter(sens_frr_plot, Age == 25), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,5.5), breaks = seq(0,5,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("C: Males and females aged 25") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_25)

plot_30 <- ggplot(data = filter(sens_frr_plot, Age == 30), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("D: Males and females aged 30") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_30)

plot_combined <- grid.arrange(plot_15, plot_20, plot_25, plot_30, ncol = 2, nrow = 2)
ggsave("output/FigureS2_2.pdf", plot = plot_combined, w=14,h=16, dpi=2400)



# RSE on FRR 0
import::from(gridExtra, grid.arrange)
sens_frr <- read_feather("output/incidences_frr_sens_rse0.feather")

sens_frr_plot <- sens_frr %>%
  select(Age, Sex, FRR, IncW, IncW.LB, IncW.UB) %>%
  mutate(FRR = FRR*100,
         IncW = IncW*100,
         IncW.LB = IncW.LB*100,
         IncW.UB = IncW.UB*100) %>%
  filter(Age %in% c(15,20,25,30))


sens_frr_plot <- sens_frr_plot %>%
  filter(Sex != "Both") %>%
  mutate(Sex = factor(case_when(
    Sex == "Men" ~ "Males",
    Sex == "Women" ~ "Females"
  ), levels = c("Males","Females")))

plot_15 <- ggplot(data = filter(sens_frr_plot, Age == 15), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("A: Males and females aged 15") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_15)

plot_20 <- ggplot(data = filter(sens_frr_plot, Age == 20), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("B: Males and females aged 20") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_20)

plot_25 <- ggplot(data = filter(sens_frr_plot, Age == 25), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("C: Males and females aged 25") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_25)

plot_30 <- ggplot(data = filter(sens_frr_plot, Age == 30), aes(x = FRR, colour = Sex)) + #, colour = Age
  #geom_line() + #ymin = IncW.LB, ymax = IncW.UB, colour = Age
  #geom_errorbar(aes(ymin = IncW.LB, ymax = IncW.UB)) +
  geom_point(aes(y = IncW)) +
  geom_smooth(aes(y = IncW), se = F) +
  geom_smooth(aes(y = IncW.LB), se = F, linetype = "dashed") +
  geom_smooth(aes(y = IncW.UB), se = F, linetype = "dashed") +
  scale_y_continuous(limits = c(0,7.1), breaks = seq(0,7,1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("D: Males and females aged 30") +
  ylab("Incidence (cases/100PY)") + xlab("FRR (%)")
print(plot_30)

plot_combined <- grid.arrange(plot_15, plot_20, plot_25, plot_30, ncol = 2, nrow = 2)
ggsave("output/FigureS2_3.pdf", plot = plot_combined, w=14,h=16, dpi=2400)
