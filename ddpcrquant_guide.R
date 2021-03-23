###############################################################################
####              INSTRUCTIONS FILE FOR DDPCRQUANT                         ####
###############################################################################

#### EXPORT DDPCR FILES FROM QUANTASOFT

#' Before starting the file export, make sure that the PC configurations 
#' for 'decimal sign' is put to a dot "."

#' Export the amplitudefiles from your ddpcr run in the Quantasoft software
#' In Quantasoft setup mode select the wells of interest and click "options". 
#' Next, click "Export Amplitude and Cluster Data". 
#' Finally, select the ddPCR folder (where your quantasoft files are located) 
#' for storage of the csv files.


#### MAKE THE DDPCRQUANT FUNCTIONS AVAILABLE

#' Make sure you have saved your ddpcrquant_function.R file in the same folder 
#' as this intructions file

#' Make the ddpcrquant functions available by opening ddpcrquant_functions.R 
#' in RStudio

#' Set the working directory to the folder containing the two R files
#' CTRL+SHIFT+H for popup window
#' 
#' Source the functions file
source(file="./ddpcrquant_functions_qx200.R")

#' the 3 functions should have appeared in the Environment tab at the right

###############################################################################
####              ddpcRquant [ggplot edition] changelog (2020 - 10 - 16)   ####
####                                                                       ####
####                           - Stevin Wilson, Ph.D. (Clemson University) ####
###############################################################################
# - Fixed installation instructions for the dpcR package.
# - Ability to specify the seed for the random number generator which is used 
#   in the analysis.
# - The output HTML file includes information about the input parameters used 
#   in the analysis.
# - Default plotting library switched to ggplot2. Ability to use the r-base 
#   plotting library via the 'plot.engine' argument in the 
#  'ddpcr.fullquant' and the 'ddpcr.threshold' functions
# - Ability to choose theme for output plots using 'plot.theme' argument in the
#   'ddpcr.fullquant' and the 'ddpcr.threshold' functions [works only when
#   plot.engine = "ggplot"] 
# - Outputs static and interactive (requires Plotly; saves as an HTML file) 
#   scatterplots with error bars to indicate the concentration across wells and 
#   samples.
# - 'ddpcr.fullquant' function returns the sample report as a tibble in 
#   addition to saving it as a csv file.


#### INSTALL PACKAGES FOR PROPER DDPCRQUANT FUNCTIONALITY

#' Install following CRAN packages
install.packages("Rcpp", type = 'source') # without this, install_github might throw an error
install.packages("tidyverse", dependencies = TRUE) # for dplyr, tibble and ggplot
install.packages("modeest", dependencies = TRUE)
install.packages("R2HTML", dependencies = TRUE)
install.packages("evd", dependencies = TRUE)
install.packages("plotly", dependencies = TRUE) # to generate interactive plot
install.packages("htmlwidgets", dependencies = TRUE) # to save interactive plot as an widget in an html file
install.packages("devtools", dependencies = TRUE) # to install the dpcR package

#' Make them available
library("modeest")
library("R2HTML")
library("evd")
library("tidyverse")
library("plotly")
library("htmlwidgets")
library("devtools")

# install dpcR from github repo
install_github("michbur/dpcR") # dpcR not available on CRAN anymore
library("dpcR")

#set seed for reproducibility
set.seed(1553)

#### ANALYSIS BY DDPCRQUANT

#' Select your working directory (folder with ddpcr files) 
#' by pressing CTRL+SHIFT+H for popup window


## FULLY AUTOMATED: ALL ASSAYS IN THE DDPCR RUN WILL BE ANALYSED

#' Run command with default settings

sample_report_tibble <- ddpcr.fullquant(work.dir = getwd(), 
                outputfile = "output_ddpcrquant",
                assay = FALSE,
                threshold.int = 0.9995,
                reps = 100, 
                blocks = 150,
                threshold.manual = FALSE,
                ci.conc.method = "wilson", 
                ci.conc = 0.95, 
                vol.mix = 20, 
                vol.sample = 2,
                plot.engine = "ggplot",
                plot.theme = "classic",
                seed.random = 1553)

#' Arguments can be changed as followed
#' 
#'   workdir = defaults to getwd(), set working directory to folder with 
#'             the amplitude.csv files and  head file
#'             Alternatively, give correct filepath (i.e. "c:/folder/file")
#'              
#'   outputfile = defaults to "output_ddpcrquant", type in other folder name
#'   
#'   assay = defaults to FALSE, numeric (cannot exceed number of assays run)
#'           if a value is set, only this assay will be run
#'           
#'   threshold.int = defaults to 0.9995, numeric, between 0-1
#'   
#'   reps = repetition of the extreme value distribution fit
#'          defaults to 100, can be changed to 1-10000
#'          more reps = more accurate = but more processing time
#'           
#'   blocks = block size, range between 100-300  
#'   
#'   threshold.manual =  defaults to FALSE, numeric
#'                       if a value is set, it will override treshold.int 
#'                       (run without extreme value theory)
#'   
#'   ci.conc.method=  all methods from dpcR are possible: wilson (default), 
#'                    agresti-coull, exact, prop.test, profile, lrt, 
#'                    asymptotic, bayes, cloglog, logit, probit
#'                    
#'   ci.conc = defaults to 0.95, numeric, between 0-1
#'   
#'   vol.mix = defaults to 20, numeric
#'   
#'   vol.sample = defaults to 2, numeric
#'
#'   for dilution factor D = vol.mix / vol.sample
#'   
#'   plot.engine = "ggplot" to use the ggplot2 package; "base" will use the plotting library from R base
#'   
#'   plot.theme can be "light", "dark", "classic" (default) or "grey" if ggplot2 is used for plotting
#'
#'   seed.random sets the input number as the seed of the random number generator


## SEMI AUTOMATED: ALL or SELECTED ASSAYS IN THE DDPCR RUN WILL BE ANALYSED

#' Run following command with default settings to explore your ddpcr run

ddpcrquant.read(work.dir=getwd())

#' Arguments can be changed as followed
#' 
#'   workdir = defaults to getwd(), so set working directory to folder with 
#'             the amplitude.csv files and  head file
#'             Alternatively, give correct filepath (i.e. "c:/folder/file")


#' 8 objects should be loaded in your Environment tab
#' Explore by running the object name
list_assays
ntc_id
sample_id


#' Run following command with the assay number of your choice (i.e. 1)
#' Assays are numbered according to the order in the ntc_id, list_assays and 
#' sample_id objects
#' Every new assay starts with a $ sign

ddpcr.fullquant(work.dir = getwd(), 
                outputfile = "output_ddpcrquant",
                assay = 1,
                threshold.int = 0.9995,
                reps = 100, 
                blocks = 150,
                threshold.manual = FALSE,
                ci.conc.method = "wilson", 
                ci.conc = 0.95, 
                vol.mix = 20, 
                vol.sample = 2,
                plot.engine = "ggplot",
                plot.theme = "classic",
                seed.random = 1553)

#' Arguments can be changed as followed
#' 
#'   workdir = defaults to getwd(), set working directory to folder with 
#'             the amplitude.csv files and  head file
#'             Alternatively, give correct filepath (i.e. "c:/folder/file")
#'              
#'   outputfile = defaults to "output_ddpcrquant", type in other folder name
#'   
#'   assay = defaults to FALSE, numeric (cannot exceed number of assays run)
#'           if a value is set, only this assay will be run
#'           
#'   threshold.int = defaults to 0.9995, numeric, between 0-1
#'   
#'   reps = repetition of the extreme value distribution fit
#'          defaults to 100, can be changed to 1-10000
#'          more reps = more accurate = but more processing time
#'           
#'   blocks = block size, range between 100-300
#'   
#'   threshold.manual =  defaults to FALSE, numeric
#'                       if a value is set, it will override treshold.int 
#'                       (run without extreme value theory)
#'   
#'   ci.conc.method=  all methods from dpcR are possible: wilson (default), 
#'                    agresti-coull, exact, prop.test, profile, lrt, asymptotic,  
#'                    bayes, cloglog, logit, probit
#'                    
#'   ci.conc = defaults to 0.95, numeric, between 0-1
#'   
#'   vol.mix = defaults to 20, numeric
#'   
#'   vol.sample = defaults to 2, numeric
#'
#'   for dilution factor D = vol.mix / vol.sample
#'   
#'   plot.engine = "ggplot" to use the ggplot2 package; "base" will use the plotting library from R base
#'   
#'   plot.theme can be "light", "dark", "classic" (default) or "grey" if ggplot2 is used for plotting
#'
#'   seed.random sets the input number as the seed of the random number generator

## SEMI AUTOMATED + TRESHOLD CONTROL

#' Run following command with default settings to explore your ddpcr run

ddpcrquant.read(work.dir=getwd())

#' Have a look at the thresholds per assay with default settings before 
#' using the ddpcr.fullquant function
#' This function needs the ddpcrquant.read function to be run before this

ddpcrquant.threshold(threshold.int = 0.9995,
                    reps = 100, 
                    blocks = 150,
                    threshold.manual = FALSE,
                    assay = 1,
                    saveplots = FALSE,
                    outputfile = FALSE,
                    ci.conc.method ="wilson", 
                    ci.conc= 0.95,
                    vol.mix = 20,
                    vol.sample = 2,
                    plot.engine = "ggplot",
                    plot.theme = "classic",
                    seed.random = 1553)

#' For every assay 2 plots should have appeared in the Plots tab on 
#' the right side of the screen: first the merged NTC and second the threshold plot
#' 
#' If the threshold does not suit the end-user, try the following:
#' 
#' Determine which assay threshold does not look good and run ddpcrquant.threshold 
#' again with selected assay number (i.e. assay = 1) and adapt the following:
#' 1) run again with adapted threshold.int (range: 0.99-0.9999995)
#' 2) run again with adapted number of blocks (range: 100-300)
#' 3) run again with a manual threshold
#' 
#' To obtain quantification results run the ddpcr.fullquant according to your 
#' changes per assay
#' 
#' Arguments can be changed as followed
#'   
#'   threshold.int = defaults to 0.9995, numeric, between 0-1
#'   
#'   reps = repetition of the extreme value distribution fit
#'          defaults to 100, can be changed to 1-10000
#'          more reps = more accurate = but more processing time
#'           
#'   blocks = block size, range between 100-300
#'   
#'   threshold.manual =  defaults to FALSE, numeric
#'                       if a value is set, it will override treshold.int 
#'                       (run without extreme value theory)
#'                       
#'   assay = defaults to FALSE, numeric (cannot exceed number of assays run)
#'           if a value is set, only this assay will be run
#'   
#'   saveplots = defaults to FALSE, if set TRUE, plots will be saved
#'   
#'   outputfile = defaults to FALSE, name can changed (i.e. "output_ddpcrquant")
#'   
#'   ci.conc.method=  all methods from dpcR are possible: wilson (default), 
#'                    agresti-coull, exact, prop.test, profile, lrt, asymptotic,  
#'                    bayes, cloglog, logit, probit
#'                    
#'   ci.conc = defaults to 0.95, numeric, between 0-1
#'   
#'   vol.mix = defaults to 20, numeric
#'   
#'   vol.sample = defaults to 2, numeric
#'
#'   for dilution factor D = vol.mix / vol.sample
#'   
#'   plot.engine = "ggplot" to use the ggplot2 package; "base" will use the plotting library from R base
#'   
#'   plot.theme can be "light", "dark", "classic" (default) or "grey" if ggplot2 is used for plotting
#'
#'   seed.random sets the input number as the seed of the random number generator

###############################################################################
####                            END OF FILE                                ####
###############################################################################
