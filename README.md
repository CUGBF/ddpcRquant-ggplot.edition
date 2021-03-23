# ddpcRquant v2

ddpcRquant is an R implementation of the algorithm described by Trypsteen et al., 2015, to calculate the threshold value for the quantification of DNA templates in a droplet digital PCR (ddPCR) experiment.

ddpcRquant v2 utilizes ggplot2 by default to generate static plots (with option available to utilize R base plotting library). Additionally, ddpcRquant v2 generates additional static and interactive plots (using Plotly), and includes updated installation procedure for dependencies.

ddpcRquant v1 is available for download at this [link](https://ddpcrquant.ugent.be/)

## Installation and Usage

Save a copy of the R scripts in this repository, to the folder containing the head file (.csv format) and the amplitude files (.csv format). Please avoid having other files in that folder. 

The installation steps for the dependencies and the usage instructions are included in the guide R script. The scripts in this repository are intended to be used only with data generated using Bio-Rad QX200.

## ddpcRquant v2 changelog (2020 - 10 - 16)

- Fixed installation instructions for the dpcR package.
- Ability to specify the seed for the random number generator which is used in the analysis.
- The output HTML file includes information about the input parameters used in the analysis.
- Default plotting library switched to ggplot2. Ability to use the r-base plotting library via the 'plot.engine' argument in the 'ddpcr.fullquant' and the 'ddpcr.threshold' functions
- Ability to choose theme for output plots using 'plot.theme' argument in the'ddpcr.fullquant' and the 'ddpcr.threshold' functions [works only when plot.engine = "ggplot"] 
- Outputs static and interactive (requires Plotly; saves as an HTML file) scatterplots with error bars to indicate the concentration across wells and samples.
- 'ddpcr.fullquant' function returns the sample report as a tibble in addition to saving it as a csv file.

## Contributing

To suggest new features or report bugs, please open an issue. Pull requests are welcome.

## Citation

If you use this tool, please cite the ddpcRquant publication :

Trypsteen, W., Vynck, M., De Neve, J., Bonczkowski, P., Kiselinova, M., Malatinkova, E., ... & De Spiegelaere, W. (2015). ddpcRquant: threshold determination for single channel droplet digital PCR experiments. Analytical and bioanalytical chemistry, 407(19), 5827-5834.

ddpcRquant v1 is available for download at this [link](https://ddpcrquant.ugent.be/)

## License
[MIT](https://choosealicense.com/licenses/mit/)
