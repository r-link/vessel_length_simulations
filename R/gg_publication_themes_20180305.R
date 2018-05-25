################################################################################
################################################################################
########	
########  ggplot graphics and color themes for plotting
########	
################################################################################
################################################################################

## Define beautiful themes for plotting         
# overall theme (containing settings for ggplot graphics display)
theme_Publication <- function(base_size = 14, base_family = "Helvetica") { # helvetica can cause troubles, but it is nice
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size = base_size, base_family = base_family)
       + theme(plot.title = element_text(face = "bold", 
                                         size = rel(1.2), hjust = 0.5), 
               text = element_text(), 
               panel.background = element_rect(colour = NA), 
               plot.background = element_rect(colour = NA), 
               panel.border = element_rect(colour = 1), 
               axis.title = element_text(face = "bold", size = 10), # size = rel(1)
               axis.title.y = element_text(angle = 90, vjust = 0), 
               axis.title.x = element_text(vjust = 1), 
               axis.text = element_text(size = 8), 
               axis.line = element_line(colour = "black"), 
               axis.ticks = element_line(), 
               panel.grid.major = element_blank(), #element_line(colour = "#f0f0f0"), 
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
               strip.background = element_blank() #element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
          ))
      
}

# theme for fill scales
scale_fill_Publication <- function(...){
      discrete_scale("fill", "Publication", 
                     manual_pal(values = c("#7fc97f", "#386cb0", "#ef3b2c", "#fdb462", "#662506", 
                                           "#a6cee3", "#fb9a99", "#984ea3", "#ffff33", "#123456", 
                                           "#888888")),
                     guide = guide_legend(order = 1, 
                                          label.theme = element_text(size = 8, angle = 0, face = "italic")), 
                     ...)
}


# legend typeface is italic b/c of species names
scale_colour_Publication <- function(...){
      discrete_scale("colour", "Publication", 
      manual_pal(values = c("#7fc97f", "#386cb0", "#ef3b2c", "#fdb462", "#662506", 
                            "#a6cee3", "#fb9a99", "#984ea3", "#ffff33", "#123456", 
                            "#888888")), 
      guide = guide_legend(order = 1, 
                           label.theme = element_text(size = 8, angle = 0, face = "italic")), 
      ...)
}

# theme for linetype scales; legend typeface is set to be plain 
scale_linetype_Publication <- function(...){
      discrete_scale("linetype", "Publication", 
      manual_pal(values = c(1, 2)), 
      guide = guide_legend(order = 2,
                           label.theme = element_text(size = 8, angle = 0, face = "plain")), 
      ...)
}

