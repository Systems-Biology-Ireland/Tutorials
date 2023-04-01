# Proteomics analysis in R

An example for the analysis of proteomics data, starting from the output of MaxQuant, can be found in `proteomics_analysis.R`. The analysis has been tested in `R=4.2`.

The analysis is based on some function I have written for proteomics analysis myself, which can be found [here](https://github.com/PhilippJunk/proteomics_analysis) on my Github.

In order to have all dependencies installed and analysis scripts downloaded, please copy the code from [`setup.R`](https://raw.githubusercontent.com/Systems-Biology-Ireland/Tutorials/main/tutorials/Tutorial_002/setup.R) and then call the copied function from R as such: `set_up("PATH/TO/DIR")`. This will create a folder at the specified location, download all required data and R files, and install all required dependencies. Depending on how up-to-date your packages are and the speed of your internet connection, this might take a couple of minutes.

To keep this tutorial working, the analysis scripts will not download the latest version, but the latest version at the time of writing this tutorial. While I do not expect big changes to the API, if you want to use the code for your own projects, I recommend checking out the latest version [here](https://github.com/PhilippJunk/proteomics_analysis).
