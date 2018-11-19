
####################################################################
### Function to store and plot D and p-values from multiple runs ###
####################################################################

# Written by Laura Dean
# Monday 19th November 2018


# The PlotABBABABA function
PlotABBABABA <- function(fasta.file, Nruns = 1) {
  # load my modified function to store the p-value as well as D value (the original function stores only the D value)
  source("C:/Users/mbzlld/Google Drive/Post Doc/my R functions/CalcD_P_value function modified from CalcD evobiR.R")
  
  # create an object in which to store information from multiple runs
  repeats <- vector(mode = "list", length = Nruns) # modify the length of this vector to alter the number of runs
  
  # write the loop
  for(y in seq_along(repeats)) {
    repeats[[y]] <- CalcD_P_value(alignment = fasta.file,
                                  sig.test = "J", ambig = "R", block.size = 100, replicate = 100)
  }
  
  # put the D and p-values into a data frame
  tab <- as.data.frame(t(sapply(repeats, function(x) c(x$D, x$p.value))))
  colnames(tab) <- c("D", "p.value")
  
  # plot a histogram of the D values
  plot_D <-
    ggplot(tab, aes(x = D)) +
    geom_histogram(bins = 1000) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 11, colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          legend.position = c(0.5,0.85),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(limits = c(-0.1,0.2), breaks = c(-0.1,0,0.1,0.2)) +
    scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
    xlab("D") +
    ylab("Variation across 1000 runs")
  
  
  # plot a histogram of the p-values
  Plot_p_value <-
    ggplot(tab, aes(x = p.value)) +
    geom_histogram(bins = 500) +
    theme(panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.text.x = element_text(size = 11, colour = "black"),
          axis.text.y = element_text(size = 11, colour = "black"),
          legend.position = c(0.5,0.85),
          legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_x_continuous(limits = c(0,1), breaks = c(0, 0.05, 0.5, 1)) +
    scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
    xlab("p value") +
    ylab("Variation across 1000 runs")
  
  # tell the function what to return
  library(patchwork)
  return(list(tab, plot_D + Plot_p_value + plot_layout(ncol = 2, nrow = 1)))
  
}

# test the new function
#PlotABBABABA(fasta.file = "C:/Users/mbzlld/Google Drive/Post Doc/RAD data analysis/ABBA BABA analysis/populations_PGDSpider_OBSE_SCAD_OBSM_BEPA/batch_4_ordered.fasta",
#             Nruns = 2)
