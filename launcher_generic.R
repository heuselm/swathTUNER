
myPath = "C:/TEMP/SwathTuner/"

if(!"shiny" %in% installed.packages())
	install.packages("shiny")

if(!"ggplot2" %in% installed.packages())
	install.packages("ggplot2")

library("shiny")
library("ggplot2")

source(paste(myPath,"R/SwathTuner.R", sep=""))

launchGUI()


