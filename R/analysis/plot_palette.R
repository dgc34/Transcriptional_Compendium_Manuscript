library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)

########## Commonly used adjustments can be saved as gg theme objects
rotate_x_labels <- theme(axis.text.x = element_text(angle = 90, hjust = 1))

theme_clean <- theme_bw() + 
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = 'white'))

#########
color_vals <- c('#F2BE84', '#C4E8F0', '#8FD0BA', 'gray50', 
                '#F2BE84', '#C4E8F0', '#8FD0BA',
                '#F2BE84', '#C4E8F0', '#8FD0BA')
names(color_vals) <- c('B', 'V', 'V/B', 'not robust', 
                       'bacterial', 'viral', 'viral/bacterial',
                       'Bacteria', 'Virus', 'viral/bacterial')
