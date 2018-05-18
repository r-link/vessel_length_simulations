################################################################################
################################################################################
########
########          Random subsamples of vessel length distributions
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Contents ---------------------------------------------------------------------
# 1. VL_subs() - random subsamples from VL_profile objects

# 1. VL_subs() -----------------------------------------------------------------
# this function returns a random subsample of a VL_profile object:
# a number of nsub (stained/empty) vessels in each cross-section

VL_subs <- function(VL_profile, nsub){
  # calculate counts of filled vessels in function of distance from cutting plane
  subs <- with(VL_profile, 
               sapply(distances,function(x){# subset of vessels present at cut x
                 sub <- left_end <= x & right_end>= x 
                 # vessels in sub that are also on the infusion point
                 infused <- left_end[sub] <= 0
                 # count infused vessels in subsample of ncond vessels present at x
                 sum(sample(infused, nsub))
               }
               )
  )
  # prepare output
  out <- c(list(subs = subs, nsub = nsub), as.list(VL_profile))
  # return output with new class VL_subs (inheriting from VL_sample)
  structure(out, class = c("VL_subs", "VL_profile"))
}
