################################################################################
################################################################################
########
########           Sampling settings for the Monte Carlo simulations 
########                    of vessel length distributions 
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Aim: Create .csv file with the settings used for the sampling of vessel length
#      distributions

# Sections:
# 1. Prerequisites     
# 2. Preparation of sampling settings

#	1. Prerequisites -------------------------------------------------------------

## clear workspace
rm(list = ls())

## load packages
# list of packages 
pkgs <-c("tidyverse", "magrittr")     

# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[, 1])]
if (length(to_install) > 0)  for (pkg in to_install) install.packages(pkg)

# load all required packages 
for (pkg in pkgs) require(pkg, character.only = T)

# function for exponential spacing of parameter settings etc.
spaces <- function(lmin, lmax, n) lmin * (lmax / lmin) ^ ((1:n - 1) / (n - 1))


#	2. Preparation of sampling settings ------------------------------------------

## define settings for Monte Carlo experiments
# 1) number of cuts
ncuts <- c(5, 10, 15, 20, 25, 30, 40, 50, 75, 100)

# 2) number of conduits in samples
nvess <- c(100, 150, 200, 300, 500, 750, 1000, 1500, 2500, 5000)

# 3) number of repetitions
runs  <- 1:100

# 4) sampling distributions
dists <- c("exp", "erlang2", "weibull", "gamma", "lnorm")     


# expand sampling grid with these settings and reorder table
settings <- expand.grid(distribution = dists, run = runs, ncuts = ncuts, nvess = nvess)%>% 
  mutate(code = factor(paste0(substr(distribution, 1, 3), "_cut", # ID for each simulation
                              formatC(ncuts, width = 3, flag = 0), "_cond", 
                              formatC(nvess, width = 3, flag = 0), "_rep", 
                              formatC(run, width = 3, flag = 0)))) %>% 
  select(distribution, run, ncuts:code) %>%                        # reorder columns
  arrange(distribution, ncuts, nvess, run) %>%                    # reorder rows
  as.tibble()                                                     # convert to tibble

# inspect table with settings
settings

# define table with parameter settings for the different distributions
parsets <- rbind(
  tibble(distribution = "exp", 
         truemean  = seq(1, 30, length.out = 100),
         par1      = 1 / truemean,
         par1_name = "rate",
         par2      = NA,
         par2_name = NA), 
  tibble(distribution = "erlang2", 
         truemean = seq(1, 30, length.out = 100), 
         par1      = 2 / truemean,
         par1_name = "rate",
         par2      = 2,
         par2_name = "shape"), 
  tibble(distribution = "weibull", 
         truemean = rep(seq(1, 30, length.out = 10), 10), 
         par1_name = "scale",
         par2 = rep(round(spaces(0.6, 6, 10), digits = 1), each = 10), 
         par2_name = "shape",
         par1      = truemean / gamma(1 + 1 / par2))[, c(1,2,6,3:5)], 
  tibble(distribution = "gamma", 
         truemean = rep(seq(1, 30, length.out = 10), 10), 
         par1_name = "rate",
         par2      = rep(round(spaces(0.6, 20, 10), digits = 1), each = 10),
         par2_name = "shape",
         par1      = par2/truemean)[, c(1,2,6,3:5)], 
  tibble(distribution = "lnorm"      , 
         truemean = rep(seq(1, 30, length.out = 10), 10),    
         par1_name = "meanlog",
         par2      = rep(round(seq(0.2, 1.2, length.out = 10), digits = 2), each = 10),
         par2_name = "sdlog",
         par1      = log(truemean) - par2^2 / 2)[, c(1,2,6,3:5)]
  )

# add columns for run, rearrange table             
parsets1 <- mutate(parsets, run = rep(1:100, 5)) %>% 
  arrange(distribution, run) %>%
  select(distribution, run, truemean:par2_name)
#inspect
parsets1

# join tables to obtain joined dataframe with all settings for the simulated distributions
settings1 <- left_join(settings, parsets)  %>% 
  arrange(distribution, ncuts, nvess, run)    %>%
  select(code, distribution, ncuts, nvess, run, truemean:par2_name)
# inspect
settings1

# export modelling settings
write.csv(settings1, "settings/monte_carlo_settings.csv", row.names = F)            
