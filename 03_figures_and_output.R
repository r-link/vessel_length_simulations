################################################################################
################################################################################
########
########      Preparation of output and figures for the Monte Carlo simulations
########                     of vessel length measurements
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Aim: Create pre-formatted summary tables and figures summarizing the results 
#      of the Monte Carlo simulations of vessel length measurements

# Sections:
# 1. Prerequisites  
#	2. Load and prepare model output
#	3. Inspection 
#	4. Table for bias & precision 
#	5. Table for overall accuracy (percent overlap)
#	6. Table for coverage 
#	7. Table for threshold on accuracy
#	8. Figures for supplementary material

#	1. Prerequisites -------------------------------------------------------------

## load packages
# list of packages 
pkgs <-c("tidyverse", "magrittr", "grid", "ggthemes", "scales")        

# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[, 1])]
if (length(to_install) > 0)  for (pkg in to_install) install.packages(pkg)

# load all required packages 
for (pkg in pkgs) require(pkg, character.only = T)

# function for equally spaced output in results tables
nice <- function(num, digits){
  # round to required number of digits
  out <- formatC(num, digits = digits, format = "f") 
  # pad with empty space if necessary
  pad <-  sapply( max(nchar(out)) - nchar(out), function(x) paste0("", rep(" ", times = x)))
  out <- paste0(pad, out)
  # return output
  return(out)
} 

# load ggplot publication themes
source("R/06_gg_themes.R")

# set figure sizes preferred by JTB (mm)
onecol     <- 90
onehalfcol <- 140
twocol     <- 190
pageheight <- 240

#	2. Load and prepare model output ---------------------------------------------

# load results of the Monte Carlo experiments
data_exh <- read_csv("output/output_exhaustive.csv", 
                     col_types = cols(par2 = col_double()))
data_sub <- read_csv("output/output_subsamples.csv", 
                     col_types = cols(par2 = col_double())) 

# merge results to a single results table, reorder, and calculate missing variables
data <- bind_rows(mutate(data_exh, Type = "Exhaustive"),
                  mutate(data_sub, Type = "Subsamples")) %>%
  filter(!is.na(code)) %>%             # remove empty first rows
  select(code, Type, everything()) %>%
  mutate(
    # properly named and ordered variables for plots
    Model = factor(model, levels =c("exp", "erlang2","gamma","weibull","lnorm","cohen","christman"), 
                   labels = c("Exponential", "Erlang(2)","Gamma","Weibull","Log-normal","Cohen","Christman"),
                   ordered =T),
    Distribution = factor(distribution, levels =c("exp", "erlang2","gamma","weibull","lnorm"), 
                          labels = c("Exponential", "Erlang(2)","Gamma","Weibull","Log-normal"),
                          ordered =T),
    SRE =    (mean_est - mean_true)/mean_true,                   # get signed relative error
    lmean_best = coalesce(lmean, lmeanquad),                     # get bet available confidence intervals for all models 
    umean_best = coalesce(umean, umeanquad),                     # (spline-based if available, quadratic approximation if not)
    included   = mean_true > lmean_best & mean_true < umean_best # get indicator for estimates that are within the CI
  ) %>%
  arrange(Type, distribution, ncuts, nvess, run, model) 
data


#	3. Inspection ----------------------------------------------------------------
# number of observations that did not converge
sum(1 - data$converged) # 6 models did not converge

# which models did not converge?
filter(data, !converged) # only the (only marginally relevant) Christman models failed

# how many re-runs with new starting values were necessary to get the models running?
table(data$t)

# how often did the calculation of the overlap coefficient fail?
sum(is.na(data$OVL))                      # 19 times
filter(data, is.na(OVL)) %$% table(model) # only for christman and weibull models


#	4. Calculate average accuracies -----------------------------------------------
# averages on distribution level 
data_dist <- data %>% 
  filter(nvess >= 500) %>% # only compare the samples that were analyzed with both methods
  group_by(Type, distribution, Model) %>%    # specify desired groups
  summarize(MSRE  = mean(SRE, na.rm = TRUE), # calculate average mean signed relative deviation
            SDSRE =   sd(SRE, na.rm = TRUE), # calculate standard deviation of MSRE (= sd of relative estimates)
            OVL   = mean(OVL, na.rm = TRUE)) # calculate average OVL

# overall averages for each model
data_mod  <- data %>% 
  filter(nvess >= 500) %>% # only compare the samples that were analyzed with both methods
  group_by(Type, Model) %>%    # specify desired groups
  summarize(MSRE  = mean(SRE, na.rm = TRUE), # calculate average mean signed relative deviation
            SDSRE =   sd(SRE, na.rm = TRUE), # calculate standard deviation of MSRE (= sd of relative estimates)
            OVL   = mean(OVL, na.rm = TRUE)) %>%  # calculate average OVL
  mutate(distribution = "All") %>%           # get values for distribution
  select(Type, distribution, everything())   # change order of columns

# combine averages to single table
data_summarized <- bind_rows(data_dist, data_mod) %>% 
  ungroup %>%
  mutate(
    # properly named and ordered variables for results tables (including average over all Distributions)
    Distribution = factor(distribution, levels =c("exp", "erlang2","gamma","weibull","lnorm", "All"), 
                          labels = c("Exponential", "Erlang(2)","Gamma","Weibull","Log-normal", "All"),
                          ordered =T)) %>%
  select(-distribution)


#	5. Table for bias & precision -------------------------------------------------
# create table for MSRE +- SDSRE
(table_bias_prec <- data_summarized %>%
   mutate(out = paste(nice(MSRE, 2), "±", nice(SDSRE, 2))) %>% # format output
   select(-MSRE, - SDSRE, - OVL) %>% # select relevant columns
   spread(key = Model, value = out)) # spread table to wide format

# export table for bias & precision
write.csv(table_bias_prec, "output/table_bias_precision.csv", row.names = FALSE)


#	5. Table for overall accuracy (percent overlap) -------------------------------
(table_OVL <- data_summarized %>%
   mutate(OVL = nice(MSRE, 3)) %>%   # format output
   select(-MSRE, - SDSRE) %>%        # select relevant columns
   spread(key = Model, value = OVL)) # spread table to wide format

# export table for bias & precision
write.csv(table_OVL, "output/table_overlap.csv", row.names = FALSE)


#	6. Table for coverage ----------------------------------------------------------
# aggregate data
(table_coverage <- data %>%
  filter(model == distribution,  # only for models fit on the correct distribution
         nvess >= 500) %>%       # ...and fit onto the same samples
  group_by(Type, Distribution) %>%  # group by distribution 
  summarize(coverage = mean(included, na.rm = TRUE)) %>% # get averages of coverage
  spread(key = Distribution, value = coverage) )  # spread table to wide format
  
# export table for coverage
write.csv(table_coverage, "output/table_coverage.csv", row.names = FALSE)


#	7. Table for threshold on accuracy --------------------------------------------
# this table is for internal use only: a table of accuracies at ncut == 10 for 
# recommendations of ncond for a specific target accuracy
table_thresh <- data %>% 
  filter(model == distribution,
         ncuts  == 10) %>%
  mutate(nvess = ifelse(Type == "Subsamples", nvess/10, nvess)) %>%
  group_by(Type, Distribution, nvess) %>%
  summarize(OVL   = mean(OVL, na.rm = TRUE),
            MSRE  = mean(SRE, na.rm = TRUE),
            SDSRE =   sd(SRE, na.rm = TRUE))

# export table for accuracy thresholds
write.csv(table_thresh, "output/accuracy_tresholds.csv", row.names = FALSE)

table_thresh %>% 
  group_by(Type, Distribution) %>%
  summarize(thresh = min(nvess[OVL >= 90]))


#	8. Figures for supplementary material -----------------------------------------
# aggregate data for plotting
plotdat <- data %>%
  filter(!(model %in% c("christman", "cohen"))) %>%      # exclude literature models from
                                                         # figures in supplementary material 
                                                         # (as requested by Reviewer 1)
  group_by(Type, ncuts, nvess, Model, Distribution) %>%  # group by relevant variables
  summarize(MSRE     = mean(SRE, na.rm = TRUE),          # summarize variables to plot
            OVL      = mean(OVL / 100, na.rm = TRUE),      
            coverage = mean(included, na.rm = TRUE)) %>%
  mutate(Ncuts = factor(paste(ncuts, "cuts"),            # create new variable for number of cuts (for facetting)
                        levels = (paste(c(5, 10, 15, 20, 25, 30, 40, 50, 75, 100),"cuts")),
                        ordered =T))

# a) bias vs ncuts:  exhaustive sampling
plotdat %>% 
  filter(Type == "Exhaustive") %>%
  ggplot(aes(x = nvess / 1000, y = MSRE, col = Model)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Mean signed relative error")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.1_bias_exh_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.1_bias_exh_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)

# b) bias vs ncuts: subsampling estimator
plotdat %>% 
  filter(Type == "Subsamples") %>% 
  ggplot(aes(x = nvess / 10000, y = MSRE, col = Model)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Mean signed relative error")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.2_bias_sub_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.2_bias_sub_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)

# c) percent overlap: exhaustive sampling
plotdat %>% 
  filter(Type == "Exhaustive") %>%
  ggplot(aes(x = nvess / 1000, y = OVL, col = Model)) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Average overlapping coefficent")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.3_overlap_exh_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.3_overlap_exh_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)

# d) percent overlap: subsampling estimator
plotdat %>% 
  filter(Type == "Subsamples") %>%
  ggplot(aes(x = nvess / 10000, y = OVL, col = Model)) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Average overlapping coefficent")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.4_overlap_sub_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.4_overlap_sub_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)

# e) coverage probabilities: exhaustive sampling
plotdat %>% 
  filter(Type == "Exhaustive") %>%
  ggplot(aes(x = nvess / 1000, y = coverage, col = Model)) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Coverage probability of confidence intervals")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.5_coverage_exh_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.5_coverage_exh_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)

# f) coverage probabilities: subsampling estimator
plotdat %>% 
  filter(Type == "Subsamples") %>%
  ggplot(aes(x = nvess / 10000, y = coverage, col = Model)) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_line(alpha = 0.8, lwd = 0.4)+
  facet_grid(Distribution ~ Ncuts,  scales = "free_y") +
  theme_Publication() + 
  scale_colour_Publication() +
  xlab("Number of vessels [10⁴]") +
  ylab("Coverage probability of confidence intervals")+  
  theme(plot.margin=unit(c(0, 0, 0, 0),"mm"),
        legend.position = "bottom",
        axis.text = element_text(size = 7),
        axis.title = element_text(vjust = -0.2, size = 11))

# export graphic
ggsave("figures/screen/figS2.6_coverage_sub_fullpage_300dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 300)
ggsave("figures/print/figS2.6_coverage_sub_fullpage_1000dpi.tiff", width = pageheight, height = twocol, units = "mm", dpi = 1000)
