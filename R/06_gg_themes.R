################################################################################
################################################################################
########
########            ggplot graphics and color themes for plotting
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# this script is based on the standard ggplot settings I use in most of my 
# projects

# Contents ---------------------------------------------------------------------
# 1. Overall theme 
# 2. Theme for colour scales 

# 1. Overall theme -------------------------------------------------------------
#  settings for ggplot graphics display
theme_Publication <- function(base_size = 14, base_family = "Helvetica") {
  (theme_foundation(base_size = base_size, base_family = base_family)
   + theme(plot.title = element_text(face = "bold", 
                                     size = rel(1.2), hjust = 0.5), 
           text = element_text(), 
           panel.background = element_rect(colour = NA), 
           plot.background = element_rect(colour = NA), 
           panel.border = element_rect(colour = 1), 
           axis.title = element_text(face = "bold", size = 10), 
           axis.title.y = element_text(angle = 90, vjust = 0), 
           axis.title.x = element_text(vjust = 1), 
           axis.text = element_text(size = 8), 
           axis.line = element_line(colour = "black"), 
           axis.ticks = element_line(), 
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), 
           legend.key = element_rect(colour = NA), 
           legend.key.height = unit(0.4, "cm"), 
           legend.title = element_text(face = "bold", size = 10), 
           legend.text  = element_text(size = 9), 
           legend.background = element_rect(fill = NA), 
           legend.spacing = unit(0, "cm"), 
           plot.margin = unit(c(1, 1, 2, 2), "mm"), 
           strip.text = element_text(face = "bold", size = 10), 
           panel.spacing = unit(2, "mm"), 
           strip.background = element_blank()
   ))
  
}

# 2. Theme for colour scales ---------------------------------------------------
scale_colour_Publication <- function(...){
  discrete_scale("colour", "Publication", 
                 manual_pal(values = c("#7fc97f", "#386cb0", "#ef3b2c", "#fdb462", "#662506", 
                                       "#a6cee3", "#fb9a99", "#984ea3", "#ffff33", "#123456", 
                                       "#888888")), 
                 guide = guide_legend(order = 1, 
                                      label.theme = element_text(size = 8, angle = 0)), 
                 ...)
}


