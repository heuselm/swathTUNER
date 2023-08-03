
INSTRUCTIONS
------------

- R must be installed in your computer.

- For Windows users: double-click the file "launcher_Windows.bat" in the folder "SwathTuner" to execute the program.
If the browser cannot display the page, refresh or close the windows (browser and command promt) and try again.
	Check the path to the R instalation in your computer and correct it:
		Open the file "launcher_Windows.bat" in a text editor and 
		update the text "C:\Program Files\R\R-3.1.1\bin\x64\".

- For any operating system: execute the R script "launcher_generic.R".
	Correct the variable "myPath" with the path to the folder "SwathTuner" in you computer:
		myPath = "D:/Program/SwathTuner/"

- If you don't have the required packages (shiny and ggplot2) a small window should open 
and you can select a server from a list, for example Switzerland.

- You can use the example input file in "SwathTuner/inst/extdata/LC-MS_peaks.tsv"
or create a new one using PeakView or any other LC-MS peak detection software 
(just change the column names accordingly, see the help (?) in the GUI).
