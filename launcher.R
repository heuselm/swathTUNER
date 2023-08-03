
if(!"shiny" %in% installed.packages())
	install.packages("shiny")

if(!"ggplot2" %in% installed.packages())
	install.packages("ggplot2")

library("shiny")
library("ggplot2")

myPath = getwd()

source("R/SwathTuner.R")

launchGUI()


