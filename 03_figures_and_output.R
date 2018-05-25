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
#	7. Figures for supplementary material

#	1. Prerequisites -------------------------------------------------------------

## load packages
# list of packages 
pkgs <-c("tidyverse", "magrittr")        

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
    SRE =    (mean_est - mean_true)/mean_true #   signed relative error
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
  mutate(lmean_best = coalesce(lmean, lmeanquad), # get bet available confidence intervals for all models 
         umean_best = coalesce(umean, umeanquad), # (spline-based if available, quadratic approximation if not)
         included   = mean_true > lmean_best & mean_true < umean_best) %>% # get indicator for estimates that are within the CI
  group_by(Type, Distribution) %>%  # group by distribution 
  summarize(coverage = mean(included, na.rm = TRUE)) %>% # get averages of coverage
  spread(key = Distribution, value = coverage) )  # spread table to wide format
  
# export table for bias & precision
write.csv(table_coverage, "output/table_coverage.csv", row.names = FALSE)


#	7. Figures for supplementary material -----------------------------------------