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

library("modeest")
library("dpcR")
library("R2HTML")
library("evd")
library("tidyverse")
library("plotly")
library("htmlwidgets")

#ddpcrquant.read()
#list_assays
#sample_id
#ddpcrquant.threshold(assay = 1, vol.temp=4, threshold.int = 0.9995)
#ddpcr.fullquant(assay =1, vol.temp=4, threshold.int = 0.9995, outputfile="RU5_9995")

#ddpcrquant.threshold(assay = 2, vol.temp=2, threshold.int = 0.9995)
#ddpcr.fullquant(assay =2, vol.temp=2, threshold.int = 0.9999, outputfile="RnaseP_9999")

ddpcr.fullquant <- function(work.dir = getwd(), 
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
                            plot.theme = "classic" ,
                            seed.random = 1553){
  ### vol.mix = volume of the mix in uL (used in the determination of concentation; dilution factor D = vol.mix / vol.sample)
  ### vol.sample = volume of the sample in uL 
  
  #set input folder
  setwd(work.dir)
  
  #set seed for reproducibility
  set.seed(seed.random)
  
  #create outputfolder
  dir.create(outputfile, showWarnings = TRUE, 
            recursive = TRUE, mode = "0777")
  output <- paste(outputfile,"/",sep="")
  
  #identify the files of the ddpcr run, make sure only ddpcr.csv files are in this folder
  ddpcrfiles <- list.files(pattern="*.csv")
  headfile <- ddpcrfiles[-grep("Amplitude", ddpcrfiles)][1]
  amplitudefiles <- list.files(pattern="*Amplitude.csv")
  
  ###### READ THE HEADFILE AND FILTER COLUMNS ######
  #read the headfile
  head_df <- read.csv(file=headfile, header = TRUE, row.names=NULL)
  #colnames(head) <- c(colnames(head)[2:dim(head)[2]], "NA") ### excluded since this line removed the column name well

  head_filtered <- head_df[,c("Well","Sample","TargetType","Target")]
  colnames(head_filtered) <- c("Well","Sample","TypeAssay","Assay")
  
  ####MAKE SUBLISTS OF THE DIFFERENT ASSAYS AND IDENTIFY THE NTC ANS THE SAMPLES######
  #look for the different assays/genes
  assays <- unique(head_filtered$Assay)  ### vector containing names of target genes
  
  #make sublist for the different assays
  list_assays=list()
  for(i in 1:length(assays)){
    name <- as.character(assays[i])
    list_assays[[name]] <- head_filtered[head_filtered$Assay==assays[i],]
  }
  
  #ntc ids per assay
  ntc=list()
  ntc_id=list()
  for(j in 1:length(list_assays)){
    list_assays_tmp <- list_assays[[j]]
    name <- as.character(assays[j])
    ntc[[name]] <- list_assays_tmp[grep("NTC",list_assays_tmp$Sample, 
                                        ignore.case = TRUE),]
    ntc_id[[name]] = ntc[[j]]$Well
  }
  
  #sample ids per assay
  sample_list=list()
  sample_id=list()
  for(j in 1:length(list_assays)){
    list_assays_tmp <- list_assays[[j]]
    name <- as.character(assays[j])
    sample_list[[name]] <- list_assays_tmp[-grep("NTC",list_assays_tmp$Sample, ignore.case = TRUE),]
    sample_id[[name]] = sample_list[[j]]$Well
  }
  
  #SET PARAMETERS
  set.seed(seed.random)
  ###reps  number of simulations
  split<-blocks #number of blocks (optimized number?)
  cutoff.quantile <- threshold.int
  
  #######START OF THE AUTOMATED ANALYSIS############
  tmp_ntc <- regexpr("Amplitude.csv", amplitudefiles) ### vector containing the position where the pattern "Amplitude.csv" occurs in each filename. helpful to extract the well label from filenames
  tmp_sample <- regexpr("Amplitude.csv", amplitudefiles)
  length_assays <- length(assays) ### number of target genes
  
  l <- assay
  
  if(l != FALSE){ length_assays <- l
  } else {l <- 1}
  
  for(l in l:length_assays){ 
    
    #Defining variables
    ntc_file <- c()
    ntc_file_corr <- c()
    channel <- c()
    report_ntc <- c()
    report_sample <- c()
    html_report <- HTMLInitFile(file.path(output),filename=paste("Report_",as.character(assays[[l]])), BackGroundColor="#FFFFFF")
    
    HTML("ddpcRquant v1.1 INPUT", file=html_report)
    HTML("Function : ddpcr.fullquant", file=html_report, append=TRUE)
    HTML(paste("Working directory :", work.dir, sep=" "), file=html_report, append=TRUE)
    HTML(paste("Name of the output folder:", outputfile, sep=" "), file=html_report, append=TRUE)
    HTML(paste("Assay :", as.character(assay), sep=" "), file=html_report, append=TRUE)
    HTML(paste("threshold.int :", threshold.int, sep=" "), file=html_report, append=TRUE)
    HTML(paste("reps :", reps, sep=" "), file=html_report, append=TRUE)
    HTML(paste("blocks :", blocks, sep=" "), file=html_report, append=TRUE)
    HTML(paste("threshold.manual :", as.character(threshold.manual), sep=" "), file=html_report, append=TRUE)
    HTML(paste("ci.conc.method :", ci.conc.method, sep=" "), file=html_report, append=TRUE)
    HTML(paste("ci.conc :", ci.conc, sep=" "), file=html_report, append=TRUE)
    HTML(paste("vol.mix :", vol.mix, sep=" "), file=html_report, append=TRUE)
    HTML(paste("vol.sample :", vol.sample, sep=" "), file=html_report, append=TRUE)
    HTML(paste("plot.engine :", plot.engine, sep=" "), file=html_report, append=TRUE)
    if (plot.engine != "base"){
      HTML(paste("plot.theme (if plot.engine = ggplot) :", plot.theme, sep=" "), file=html_report, append=TRUE)
    }
    
    HTML("                                ", file=html_report, append=TRUE)
    HTML("********************************", file=html_report, append=TRUE)
    HTML("--------------------------------", file=html_report, append=TRUE)
    HTML("********************************", file=html_report, append=TRUE)
    HTML("                                ", file=html_report, append=TRUE)
    
    
    HTML("ddpcRquant v1.1 OUTPUT", file=html_report, append=TRUE)
    HTML(paste("OUTPUT FOR ASSAY:", assays[l], sep=" "), file=html_report, append=TRUE)
    
    ###if(unique(list_assays[[l]]$TypeAssay)=="Ch1Unknown"){
    ###  channel <- 1
    ###} else{
    ###  channel <- 2
    ###}
    ### replaced with a more robust option
    ### if "Ch1Unknown" is in the column "TypeAssay", then consider values in column 1; otherwise column 2

    if ("Ch1Unknown" %in% unique(list_assays[[l]]$TypeAssay)){
      channel <- 1
    } else{
      channel <- 2      
    }

    ####NTC PROCESSING####
    for(k in 1:length(ntc_id[[l]])){
      for(j in 1:length(amplitudefiles)){
        
        match <- substr(amplitudefiles[j], tmp_ntc-4, tmp_ntc-2)  ### extract the well labels from the filenames of amplitude files
        
        if (match == ntc_id[[l]][k]){
          ### if the well label in the amplitude filename matches that of a NTC well label
          #######READ NTCS#######
          ntc_file_tmp <- read.csv(file=amplitudefiles[j], header = T)
          
          #######BASELINE TO ZERO BASED ON ROBUST MODE#######
          corr_factor <- hsm(ntc_file_tmp[,channel])    ### half sample mode estimator
          ntc_file_corr <- ntc_file_tmp[,channel] - corr_factor   ### vector of length equal to the number of droplets, with values equal to the number in the column 1 or 2  minus corr_factor
          
          #######EVALUATE BASELINE#######
          #par(mfrow=c(2,2),oma=c(0,0,2,0))
          
          #plot(density(ntc_file_tmp[,channel]), main="")
          #abline(v=corr_factor, col="red")
          #title("not_baselined_density", round(corr_factor, 2))
          
          #plot(sample(ntc_file_tmp[,channel]),cex = 0.8, pch=20, 
          #     xlab=c("Droplets (#)"), ylab=c("Fluorescence"))
          #abline(h=corr_factor, col="red", )
          #title("not_baselined_dotplot", round(corr_factor, 2))
          
          #plot(density(ntc_file_corr), main="")
          #abline(v=0, col="red")
          #title("baselined_density")
          name_base <- paste("Baseline_",assays[l],"_ntc_",match, sep="")
          name_dens <- paste(output,"Baseline_", "ntc_",match,"_",assays[l],".png", sep="")
          
          if (plot.engine == "base"){
            plot(sample(ntc_file_corr),cex = 0.8, pch=20, 
                xlab=c("Droplets (#)"), ylab=c("Fluorescence"))    ### sample() does random reordeing of the vector ntc_file_corr
            abline(h=0, col="red")
            title(name_base)
            ### output filename for the plot
            dev.copy(png, file = name_dens)
            dev.off()
          } else {
            ntc_plot <- ggplot(mapping = aes(x = seq(length(ntc_file_corr)), y = sample(ntc_file_corr))) + geom_point(alpha = 0.3, color = "#56B4E9") + geom_hline(yintercept = 0, linetype = "dashed",color = "#D55E00") + 
              labs(x = "Droplets (#)", y = "Fluorescence", title = name_base)
            if (plot.theme == "light"){
              ntc_plot + theme_light()
            } else if (plot.theme == "dark"){
              ntc_plot + theme_dark()
            } else if (plot.theme == "classic"){
              ntc_plot + theme_classic()
            } else {
              ntc_plot + theme_grey()
            }
            ggsave(name_dens, width = 7, height = 7)
          }
          
          ###HTML OUTPUT###
          #HTML("BASELINING NTCS", file=html_report, append=TRUE)
          #HTMLInsertGraph(GraphFileName=paste("Baseline_", "ntc_",match,"_",assays[l],".png", sep=""), Caption="", GraphBorder=1, Align="center", file=html_report, append=TRUE)
          
          #######MERGE NTCS#######
          ntc_file <- c(ntc_file, ntc_file_corr) ### update ntc_file by concatenating ntc_file_corr
        }
      }}       
    
    ##### PLOT MERGED NTC ####

    #set seed for reproducibility
    set.seed(seed.random)
    
    random_ntc_file <- sample(ntc_file)
    density_max_merge <- 0 ### not sure what was the need of this line 
    density_max_merge <- hsm(ntc_file)
    
    #par(mfrow=c(1,2), oma=c(0,0,2,0))
    #plot(density(random_ntc_file), main="")
    #abline(v=density_max_merge, col="red")
    #title("density", density_max_merge)
    name_merged <- paste("Merged_NTC_", assays[l], "_#ntc_",length(ntc_id[[l]]), sep="")
    name_merge <- paste(output,"Merged_NTC_",assays[l],".png", sep="") ### output filename for the merged NTC plot
    
    if (plot.engine == "base"){ 
      plot(random_ntc_file,cex = 0.8, pch=20, xlab=c("Droplets (#)"), ylab=c("Fluorescence"))
      abline(h=density_max_merge, col="red" )
      #title("dotplot", round(density_max_merge))
      title(name_merged, paste("baseline = ", round(density_max_merge), sep=""))
      dev.copy(png, file = name_merge)
      dev.off()
    } else {
      merged_ntc_plot <- ggplot(mapping = aes(x = seq(length(random_ntc_file)), y = random_ntc_file)) + geom_point(alpha = 0.3, color = "#56B4E9") + geom_hline(yintercept = density_max_merge, linetype = "dashed",color = "#D55E00") + 
        labs(x = "Droplets (#)", y = "Fluorescence", title = name_merged, caption = paste("baseline = ", round(density_max_merge), sep=""))
      if (plot.theme == "light"){
        merged_ntc_plot + theme_light()
      } else if (plot.theme == "dark"){
        merged_ntc_plot + theme_dark()
      } else if (plot.theme == "classic"){
        merged_ntc_plot + theme_classic()
      } else {
        merged_ntc_plot + theme_grey()
      }
      ggsave(name_merge, width = 7, height = 7)
    }
    
    #HTML OUTPUT
    HTML("MERGED NTC PLOT", file=html_report, append=TRUE)
    HTMLInsertGraph(GraphFileName=paste("Merged_NTC_",assays[l],".png", sep=""), Caption="", GraphBorder=1, Align="center", file=html_report, append=TRUE)
    
    ##### CALCULATE THRESHOLD ####
    
    if(threshold.manual == FALSE){
  
      #set seed for reproducibility
      set.seed(seed.random)
      
      #initialize vectors
      quantgev<-array(0,reps)  ### array of 0s of length equal to reps; used to store output value of each rep
      
      #start repeated calculations
      for(k in 1:reps){
        print(k) #TRACK PROGRESS
        
        ### block creation
        random_ntc_file <- sample(random_ntc_file)
        size<-floor(length(random_ntc_file)/(split*length(ntc_id[[l]]))) #block size
        ssize<- c(rep(size,split*length(ntc_id[[l]])-(length(random_ntc_file)-size*split*length(ntc_id[[l]]))),
                  rep(size+1,length(random_ntc_file)-(size*split*length(ntc_id[[l]])))) #split remaining droplets over samples so no droplets get discarded ### split = blocks, which is an input parameter
        signal.maxima<-array(0,split) #to store the block maximums
        #droplets<-sample(droplets) #random subsampling, don't forget to set seed if reproducibility needed
        #-length(random_ntc_file)%%size
        
        ### determine maxima of subsamples
        droplets2<-sample(random_ntc_file)
        for(i in 1:split){
          signal.maxima[i]<-max(droplets2[1:ssize[i]])
          droplets2<-droplets2[-c(1:ssize[i])]
        }
        
        signal.maxima<-as.numeric(signal.maxima)
        
        
        ## fit the GEV model using ML
        droplet.fit <- try(fgev(signal.maxima), silent=TRUE) ### maximum likelihood fitting of generalized extreme value distribution
        
        #fgev(signal.maxima, start=data.frame(loc=400,scale=100,shape=0))
        
        #droplet.fit<-fgev(signal.maxima)
        quantgev[k]<-try(qgev(cutoff.quantile,droplet.fit$estimate[1],
                              droplet.fit$estimate[2],droplet.fit$estimate[3]), silent =TRUE) ### cutoff.quantile equals threshold.int, which is an input parameter
        if(quantgev[k] > max(random_ntc_file)+3000){quantgev[k] <- NA}  ### if the threshold value from a rep exceeds the maximum of random_ntc_file vector + 3000, then consider it NA
        
      }
      
      #calculate final threshold
      threshold<-mean(as.numeric(quantgev), na.rm=TRUE)  ### find mean of values in the quantgev vector while ignoring NA entries
    } else { threshold <- threshold.manual}
    
    ######POSITIVE AND NEGATIVE DROPLETS########
    pos_drops <- random_ntc_file[random_ntc_file>threshold]  ### logical vector with TRUE for droplets with value > threshold, and FALSE otherwise
    neg_drops <- random_ntc_file[random_ntc_file<=threshold]  ### logical vector with TRUE for droplets with value <= threshold, and FALSE otherwise
    nr_pos_drops <-length(pos_drops)  ### number of positive droplets
    nr_neg_drops <-length(neg_drops)  ### number of negative droplets
    tot <- length(random_ntc_file)
    
    ########PLOT#########
    ymax <- max(c(threshold,max(random_ntc_file))) + 1000
    report <- paste("pos: ",nr_pos_drops," neg: ",nr_neg_drops," tot: ",tot,"threshold: ",round(threshold,digits=2)) ### indicate the +ve and -ve numbers below the plot
    name_thres <- paste("merged_NTC_threshold_",assays[l], sep="")
    name_plot <- paste(output,"merged_NTC_threshold_",assays[l],".png", sep="")
    
    if (plot.engine == "base"){    
      plot(random_ntc_file, main="",ylim=c(min(random_ntc_file),ymax),cex = 0.8, 
          pch=20, ylab=c("Fluorescence"))
      abline(h = threshold, col=2, lty=4)
      title(name_thres,report)
      dev.copy(png, file = name_plot)
      dev.off()
    } else{
      merged_ntc_threshold_plot <- ggplot(mapping = aes(x = seq(length(random_ntc_file)), y = random_ntc_file)) + geom_point(alpha = 0.3, color = "#56B4E9") + geom_hline(yintercept = threshold, linetype = "dashed",color = "#D55E00") + 
        labs(x = "Droplets (#)", y = "Fluorescence", title = name_thres, caption = report) + ylim(min(random_ntc_file),ymax)
      if (plot.theme == "light"){
        merged_ntc_threshold_plot + theme_light()
      } else if (plot.theme == "dark"){
        merged_ntc_threshold_plot + theme_dark()
      } else if (plot.theme == "classic"){
        merged_ntc_threshold_plot + theme_classic()
      } else {
        merged_ntc_threshold_plot + theme_grey()
      }
      ggsave(name_plot, width = 7, height = 7)      
    }
    
    HTML("NTC THRESHOLD PLOT", file=html_report, append=TRUE)
    HTMLInsertGraph(GraphFileName= paste("merged_NTC_threshold_",assays[l],".png", sep=""), Caption="", GraphBorder=1, Align="center", file=html_report, append=TRUE)
    HTML("SAMPLE PLOTS", file=html_report, append=TRUE)
    
    ###PINHEIRO CONCENTRATION CALCS        
    if(nr_pos_drops != 0){
      dpcr <- dpcr_density(nr_pos_drops,tot, average = TRUE, 
                          methods = ci.conc.method, 
                          conf.level = ci.conc, plot = FALSE)
      
      conc <- (1000/0.91 * dpcr$lambda)*(vol.mix/vol.sample)   ### see Trypsteen et al., 2015; conc is the concentration in the number of copies per uL stock ; equation to find C; Vd = 0.91uL; dilution factor D = volume of mix(uL) / volume of sample(uL)
      conc_lower <- (1000/0.91 * dpcr$lower)*(vol.mix/vol.sample)
      conc_upper <- (1000/0.91 * dpcr$upper)*(vol.mix/vol.sample)
      
    } else{ 
      conc_lower <- 0
      conc_upper <- 0
      conc <- 0
    }
    
    report_ntc <- c("merged", as.character(assays[l]), 
                    "merged_NTC", "ntc", nr_pos_drops, 
                    nr_neg_drops, tot, round(conc, digits=4), 
                    round(conc_lower, digits=4), round(conc_upper, digits=4)) ### this is used as a row in the table within the HTML file
    
    
    ####SAMPLES PROCESSING####
    for(a in 1:length(sample_id[[l]])){
      for(b in 1:length(amplitudefiles)){
        
        match_sample <- substr(amplitudefiles[b], tmp_sample-4, tmp_sample-2)
        
        if (match_sample == sample_id[[l]][a]){
          ### if sample well label matches the well label in the name of an amplitude filename
          #######READ SAMPLE FILE#######
          sample_file<-read.csv(file=amplitudefiles[b], header = T)
          
          #######RANDOMIZE SAMPLE#######
          random_sample_file <- sample(sample_file[,channel])
          
          ##baseline to zero depending robust mode
          window <- subset(random_sample_file, 
                          random_sample_file < threshold + corr_factor) ### filter rows to have only the ones with values in the channel column less than the sum of threshold and corr_factor
          window.mode <- hsm(window)
          if(is.na(window.mode)){window.mode = 0}
          #new.cutoff <- window.mode + threshold
          #neg.droplets <- sum(random_sample_file<new.cutoff)
          #pos.droplets <- sum(random_sample_file>=new.cutoff)
          #random_sample_file_corr <- random_sample_file - hsm(random_sample_file)
          random_sample_file_corr <- random_sample_file - window.mode
          
          ######RETURN POSITIVE AND NEGATIVE DROPLETS########
          pos_drops <- sum(random_sample_file_corr>threshold) ### TRUE = 1, FALSE = 0; this line gives the number of droplets with values more than the threshold
          neg_drops <- sum(random_sample_file_corr<=threshold) ### TRUE = 1, FALSE = 0; this line gives the number of droplets with values less than or equal to the threshold
          
          tot <- length(random_sample_file_corr)
          
          ###PINHEIRO CONCENTRATION CALCS        
          if(pos_drops != 0){
            dpcr <- dpcr_density(pos_drops,tot, average = TRUE, 
                                methods = ci.conc.method, 
                                conf.level = ci.conc, plot = FALSE)
            conc <- (1000/0.91 * dpcr$lambda)*(vol.mix/vol.sample)
            conc_lower <- (1000/0.91 * dpcr$lower)*(vol.mix/vol.sample)
            conc_upper <- (1000/0.91 * dpcr$upper)*(vol.mix/vol.sample)
            
          } else{ 
            conc <- 0
            conc_lower <- 0
            conc_upper <- 0
          }
          
          report <- paste("pos: ",pos_drops," neg: ",neg_drops," tot: ",
                          tot," conc: ",round(conc, digits=4),"copies/uL stock") ### indicate the number of +ve, -ve droplets below the plot
          tmp_rep <- c(match_sample, as.character(assays[l]), 
                      as.character(list_assays[[l]][list_assays[[l]]$Well == match_sample ,]$Sample), "sample", 
                      pos_drops, neg_drops, tot, round(conc, digits=4), 
                      round(conc_lower, digits=4), round(conc_upper, digits=4)) ### used as row in the table within the HTML file
          
          
          #####PLOT SAMPLE#####
          name_sample <- paste(assays[l],match_sample,
                              list_assays[[l]][list_assays[[l]]$Well == match_sample,]$Sample,sep="_")
          sample_plot_name <- paste(output,"treshold_plot_",assays[l],"_",match_sample,"_",".png", sep="")
          
          if (plot.engine == "base"){
            plot(random_sample_file_corr, main="",cex = 0.8, pch=20, 
                ylim=c(min(random_sample_file_corr), max(random_sample_file_corr+1000)), 
                ylab=c("Fluorescence"))
            abline(h = threshold, col=2, lty=4)
            text(1000,threshold,round(threshold, digits=2), col=2)
            title(name_sample,report)
            dev.copy(png, file = sample_plot_name)
            dev.off()
          } else{
            sample_threshold_plot <- ggplot(mapping = aes(x = seq(length(random_sample_file_corr)), y = random_sample_file_corr)) + geom_point(alpha = 0.3, color = "#56B4E9") + geom_hline(yintercept = threshold, linetype = "dashed",color = "#D55E00") + 
              labs(x = "Droplets (#)", y = "Fluorescence", title = name_sample, caption = report) + ylim(min(random_sample_file_corr), max(random_sample_file_corr+1000)) + geom_text(aes(1000,threshold,label = round(threshold,3), vjust = -1), color = "#D55E00", fontface = "bold")
            if (plot.theme == "light"){
              sample_threshold_plot + theme_light()
            } else if (plot.theme == "dark"){
              sample_threshold_plot + theme_dark()
            } else if (plot.theme == "classic"){
              sample_threshold_plot + theme_classic()
            } else {
              sample_threshold_plot + theme_grey()
            }
            ggsave(sample_plot_name, width = 7, height = 7)            
          }
          
          #HTMLInsertGraph(GraphFileName=paste("treshold_plot_",assays[l],"_",match_sample,"_",".png", sep=""), Caption="", GraphBorder=1, Align="center", file=html_report, append=TRUE)
          
        }
        
      }
      
      report_sample <- rbind(report_sample, tmp_rep)
      
    }
    
    final_report <- rbind(report_ntc, report_sample)
    colnames(final_report) <- c("Well", "assay", "name", "type", "positive droplets", "negative droplets", "total droplets", "concentration", "lowerCI", "upperCI")
    rownames(final_report) <- c(1:(length(list_assays[[l]]$Well)-length(ntc_id[[l]])+1))
    
    sample_report_df <-  as.data.frame(report_sample)  # the function will return this dataframe as a tibble
    colnames(sample_report_df) <- c("Well", "assay", "name", "type", "positive_droplets", "negative_droplets", "total_droplets", "concentration", "lowerCI", "upperCI")
    sample_report_df[c(1,2,3,4)] <- lapply(sample_report_df[c(1,2,3,4)], as.factor)
    sample_report_df[c(5,6,7)] <- lapply(sample_report_df[c(5,6,7)], as.integer)
    sample_report_df[c(8,9,10)] <- lapply(sample_report_df[c(8,9,10)], as.numeric)
    sample_report_df <- as_tibble(sample_report_df)
    
    ### scatterplot comparing concentration across wells with error bars representing the confidence interval  
    colourCount <- length(unique(sample_report_df$name))
    samples_concentration_plot <- ggplot(data = sample_report_df, mapping = aes(x = Well, y = concentration, text = paste(paste("lower limit of the ", ci.conc * 100,"% CI : ", sep = ""), lowerCI,
                                                                                                                          paste("</br>upper limit of the ", ci.conc * 100,"% CI : ", sep = ""), upperCI,
                                                                                                                          "</br>positive droplets :", positive_droplets,
                                                                                                                          "</br>negative droplets : ", negative_droplets,
                                                                                                                          "</br>total droplets : ", total_droplets))) +
      labs(x = "well label", y = "concentration (in copies per uL of stock)", title = paste("Concentration vs Well Labels - ", assays[l], sep = ""), caption = paste("Error bars indicate ", as.integer(ci.conc * 100), "% confidence interval", sep = ""))

    if (plot.theme == "classic" || plot.engine != "ggplot"){
      samples_concentration_plot <- samples_concentration_plot + geom_point() + 
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI)) + 
        theme_classic() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))        
    } else if (plot.theme == "light"){
      samples_concentration_plot <- samples_concentration_plot + geom_point() + 
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI)) + 
        theme_light() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
    } else if (plot.theme == "dark"){
      samples_concentration_plot <- samples_concentration_plot + geom_point(color = "#56B4E9") + 
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), color = "#56B4E9") + 
        theme_dark() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      samples_concentration_plot <- samples_concentration_plot + geom_point() + 
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI)) + 
        theme_grey() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
    }
    
    samples_concentration_plot
    ggsave(paste(output,"concentration_across_wells_", assays[l], ".png", sep = ""), width = 7, height = 7)
    
    ### generative an interactive plot using Plotly and save it as an html file
    samples_concentration_plotly <- ggplotly(samples_concentration_plot , tooltip = c("x", "y", "text"))
    samples_concentration_plotly <- samples_concentration_plotly %>% 
      layout(annotations = 
               list(x = 1, y = -0.1, text = paste("Error bars indicate ", as.integer(ci.conc * 100), "% confidence interval.", sep = ""), 
                    showarrow = F, xref='paper', yref='paper', 
                    xanchor='right', yanchor='auto', xshift=0, yshift=0,
                    font=list(size=15, color="#56B4E9"))
      )
    htmlwidgets::saveWidget(samples_concentration_plotly, file = file.path(getwd(), paste(output,"concentration_across_wells_", assays[l], ".html", sep = "")))

    ### Group replicates of different samples and plot the concentration across the samples. error bars represent std deviation    
    final_conc_plot <- sample_report_df %>%
      group_by(name) %>%
      summarise(mean_conc = mean(concentration), stddev_conc = sd(concentration), n_replicates = n()) %>%
      ggplot(aes(x = name, y = mean_conc, text = paste("standard deviation :", stddev_conc,
                                                      "</br>number of replicates :", n_replicates)))  +
      labs(x = "sample name", y = "mean concentration (in copies per uL of stock)", title = paste("Concentration across samples - ", assays[l], sep = ""), caption = "Error bars represent the standard deviation")
    
    if (plot.theme == "classic" || plot.engine != "ggplot"){
      final_conc_plot <- final_conc_plot + geom_point() + 
        geom_errorbar(aes(ymin = mean_conc - stddev_conc, ymax = mean_conc + stddev_conc)) + 
        theme_classic()    
    } else if (plot.theme == "light"){
      final_conc_plot <- final_conc_plot + geom_point() + 
        geom_errorbar(aes(ymin = mean_conc - stddev_conc, ymax = mean_conc + stddev_conc)) + 
        theme_light()
    } else if (plot.theme == "dark"){
      final_conc_plot <- final_conc_plot + geom_point(color = "#56B4E9") + 
        geom_errorbar(aes(ymin = mean_conc - stddev_conc, ymax = mean_conc + stddev_conc), color = "#56B4E9") + 
        theme_dark()
    } else {
      final_conc_plot <- final_conc_plot + geom_point() + 
        geom_errorbar(aes(ymin = mean_conc - stddev_conc, ymax = mean_conc + stddev_conc)) + 
        theme_grey()
    }
    
    final_conc_plot
    ggsave(paste(output,"concentration_across_samples_", assays[l], ".png", sep = ""), width = 7, height = 7)

    ### generative an interactive plot using Plotly and save it as an html file
    final_conc_plotly <- ggplotly(final_conc_plot , tooltip = c("x", "y", "text"))
    final_conc_plotly <- final_conc_plotly %>% 
      layout(annotations = 
                    list(x = 1, y = -0.1, text = "Error bars represent the standard deviation.", 
                    showarrow = F, xref='paper', yref='paper', 
                    xanchor='right', yanchor='auto', xshift=0, yshift=0,
                    font=list(size=15, color= "#56B4E9"))
      )
    htmlwidgets::saveWidget(final_conc_plotly, file = file.path(getwd(), paste(output,"concentration_across_samples_", assays[l], ".html", sep = "")))  
    
    write.csv(final_report, file=paste(output, "final_report_", assays[l], ".csv", sep=""))
    
    HTML("FINAL REPORT", file=html_report, append=TRUE)
    HTML(final_report, file=html_report, append=TRUE)
    
    browseURL(paste(file.path(output, paste("Report_",assays[l], ".html", sep=""))))
    
  }
  return(sample_report_df)
}


#ddpcrquant.read()

ddpcrquant.read <- function(work.dir=getwd()){
  
  #set input folder
  setwd(work.dir)
  
  #identify the files of the ddpcr run, make sure only ddpcr.csv files are in this folder
  ddpcrfiles <<- list.files(pattern="*.csv")  ### <<- means non-local assignment
  headfile <<- ddpcrfiles[-grep("Amplitude", ddpcrfiles)][1]
  amplitudefiles <<- list.files(pattern="*Amplitude.csv")
  
  ###### READ THE HEADFILE AND FILTER COLUMNS ######
  #read the headfile
  head_df <- read.csv(file=headfile, header = TRUE, row.names=NULL)
  #colnames(head) <- c(colnames(head)[2:dim(head)[2]], "NA") ### excluded since this line removed the column name well
  
  head_filtered <- head_df[,c("Well","Sample","TargetType","Target")]
  colnames(head_filtered) <- c("Well","Sample","TypeAssay","Assay")
  
  
  ####MAKE SUBLISTS OF THE DIFFERENT ASSAYS AND IDENTIFY THE NTC ANS THE SAMPLES######
  #look for the different assays/genes
  assays <<- unique(head_filtered$Assay)
  
  #make sublist for the different assays
  list_assays=list()
  for(i in 1:length(assays)){
    name <- as.character(assays[i])
    list_assays[[name]] <- head_filtered[head_filtered$Assay==assays[i],]
  }
  list_assays <<- list_assays
  
  #ntc ids per assay
  ntc=list()
  ntc_id=list()
  for(j in 1:length(list_assays)){
    list_assays_tmp <- list_assays[[j]]
    name <- as.character(assays[j])
    ntc[[name]] <- list_assays_tmp[grep("NTC",list_assays_tmp$Sample, 
                                        ignore.case = TRUE),]
    ntc_id[[name]] = ntc[[j]]$Well
  }
  
  ntc_id <<- ntc_id
  
  #sample ids per assay
  sample_list=list()
  sample_id=list()
  for(j in 1:length(list_assays)){
    list_assays_tmp <- list_assays[[j]]
    name <- as.character(assays[j])
    sample_list[[name]] <- list_assays_tmp[-grep("NTC",list_assays_tmp$Sample, ignore.case = TRUE),]
    sample_id[[name]] = sample_list[[j]]$Well
  }
  sample_id <<- sample_id
  rm(i, name)  ### remove i and name objects
}

#assays
#ddpcrquant.threshold()

ddpcrquant.threshold <- function(threshold.int = 0.9995,
                                reps = 100, 
                                blocks = 150,
                                threshold.manual = FALSE,
                                assay = 1,
                                saveplots = FALSE,
                                outputfile = FALSE,
                                ci.conc.method ="wilson", 
                                ci.conc= 0.95,
                                vol.mix = 20,
                                vol.sample = 4 ,
                                plot.engine = "ggplot",
                                plot.theme = "classic",
                                seed.random = 1553){
  
  if(saveplots == TRUE){
    #create outputfolder
    dir.create(outputfile, showWarnings = TRUE, 
              recursive = TRUE, mode = "0777")
    output <- paste(outputfile,"/",sep="")
  }
  
  #SET PARAMETERS
  set.seed(seed.random)
  #reps <- 10 #100 number of simulations
  split <- blocks #number of blocks
  cutoff.quantile = threshold.int
  
  #######START OF THE AUTOMATED ANALYSIS############
  tmp_ntc <- regexpr("Amplitude.csv", amplitudefiles)
  report_thresh <- c()
  
  length_assays <- length(assays)

  l <- assay
  
  if(l != FALSE){ length_assays <- l
  } else {l <- 1}
  
  for(l in l:length_assays){

    #Defining variables
    ntc_file <- c()
    ntc_file_corr <- c()
    channel <- c()
    report_ntc <- c()
    
    
    ###if(unique(list_assays[[l]]$TypeAssay)=="Ch1Unknown"){
    ###  channel <- 1
    ###} else{
    ###  channel <- 2
    ###}
    ### replaced with a more robust option
    ### if "Ch1Unknown" is in the column "TypeAssay", then consider values in column 1; otherwise column 2
    
    if ("Ch1Unknown" %in% unique(list_assays[[l]]$TypeAssay)){
      channel <- 1
    } else{
      channel <- 2      
    }
    
    ####NTC PROCESSING####
    for(k in 1:length(ntc_id[[l]])){
      for(j in 1:length(amplitudefiles)){
        
        match <- substr(amplitudefiles[j], tmp_ntc-4, tmp_ntc-2)
        
        if (match == ntc_id[[l]][k]){
  
          #######READ NTCS#######
          ntc_file_tmp <- read.csv(file=amplitudefiles[j], header = T)
          
          #######BASELINE TO ZERO BASED ON ROBUST MODE#######
          corr_factor <- hsm(ntc_file_tmp[,channel])
          ntc_file_corr <- ntc_file_tmp[,channel] - corr_factor
          
          #######MERGE NTCS#######
          ntc_file <- c(ntc_file, ntc_file_corr)
        }
      }}       
    
    ##### PLOT MERGED NTC ####

    #set seed for reproducibility
    set.seed(seed.random)
    
    random_ntc_file <<- sample(ntc_file)
    density_max_merge <- 0
    density_max_merge <- hsm(ntc_file)
    
    #plot(random_ntc_file,cex = 0.8, pch=20, xlab=c("Droplets (#)"), ylab=c("Fluorescence"))
    #abline(h=density_max_merge, col="red", )
    #name_merged <- paste("Merged_NTC_", assays[l], "_#ntc_",length(ntc_id[[l]]), sep="")
    #title(name_merged, paste("baseline = ", round(density_max_merge), sep=""))
    
    #if(saveplots == TRUE){
    #  name_merge <- paste(output,"Merged_NTC_",assays[l],".png", sep="")
    #  dev.copy(png, file = name_merge)
    #  dev.off()
    #}
    ##### CALCULATE THRESHOLD ####
    
    if(threshold.manual == FALSE){
      
      #initialize vectors
      quantgev<-array(0,reps)
      
      #start repeated calculations
      for(k in 1:reps){
        print(k) #TRACK PROGRESS
        
        ### block creation
        random_ntc_file <- sample(random_ntc_file)
        size<-floor(length(random_ntc_file)/(split*length(ntc_id[[l]]))) #block size
        ssize<- c(rep(size,split*length(ntc_id[[l]])-(length(random_ntc_file)-size*split*length(ntc_id[[l]]))),
                  rep(size+1,length(random_ntc_file)-(size*split*length(ntc_id[[l]])))) #split remaining droplets over samples so no droplets get discarded
        signal.maxima<-array(0,split) #to store the block maximums
        
        ### determine maxima of subsamples
        droplets2<-sample(random_ntc_file)
        for(i in 1:split){
          signal.maxima[i]<-max(droplets2[1:ssize[i]])
          droplets2<-droplets2[-c(1:ssize[i])]
        }
        
        signal.maxima<-as.numeric(signal.maxima)
        
        
        ## fit the GEV model using ML
        droplet.fit <- try(fgev(signal.maxima), silent=TRUE)
        
        quantgev[k]<-try(qgev(cutoff.quantile,droplet.fit$estimate[1],
                              droplet.fit$estimate[2],droplet.fit$estimate[3]), silent =TRUE)
        if(quantgev[k] > max(random_ntc_file)+3000){quantgev[k] <- NA}
        
      }
      
      #calculate final threshold
      threshold<<-mean(as.numeric(quantgev), na.rm=TRUE)
    } else { threshold <- threshold.manual}
    

    ######POSITIVE AND NEGATIVE DROPLETS########
    pos_drops <- random_ntc_file[random_ntc_file>threshold]
    neg_drops <- random_ntc_file[random_ntc_file<=threshold]
    nr_pos_drops <-length(pos_drops)
    nr_neg_drops <-length(neg_drops)
    tot <- length(random_ntc_file)
    
    ########PLOT#########
    ymax <<- max(c(threshold,max(random_ntc_file))) + 1000
    report <- paste("pos: ",nr_pos_drops," neg: ",nr_neg_drops," tot: ",tot,"threshold: ",round(threshold,digits=2))
    name_thres <- paste("merged_ntc_threshold_",assays[l], sep="")
    name_plot <- paste(output,"merged_treshold",assays[l],".png", sep="")
    
    if (plot.engine == "base"){
      plot(random_ntc_file, main="",ylim=c(min(random_ntc_file),ymax),cex = 0.8, pch=20)
      abline(h = threshold, col=2, lty=4)
      title(name_thres,report)
      if(saveplots == TRUE){
        dev.copy(png, file = name_plot)
        dev.off()
      }
    } else {
      manual_merged_ntc_threshold_plot <- ggplot(mapping = aes(x = seq(length(random_ntc_file)), y = random_ntc_file)) + geom_point(alpha = 0.3, color = "#56B4E9") + geom_hline(yintercept = threshold, linetype = "dashed",color = "#D55E00") + 
        labs(x = "Droplets (#)", y = "Fluorescence", title = name_thres, caption = report) + ylim(min(random_ntc_file),ymax)
      if (plot.theme == "light"){
        manual_merged_ntc_threshold_plot + theme_light()
      } else if (plot.theme == "dark"){
        manual_merged_ntc_threshold_plot + theme_dark()
      } else if (plot.theme == "classic"){
        manual_merged_ntc_threshold_plot + theme_classic()
      } else {
        manual_merged_ntc_threshold_plot + theme_grey()
      }
      ggsave(name_plot, width = 7, height = 7) 
    }
    
    ###PINHEIRO CONCENTRATION CALCS        
    if(nr_pos_drops != 0){
      dpcr <- dpcr_density(nr_pos_drops,tot, average = TRUE, 
                          methods = ci.conc.method, 
                          conf.level = ci.conc, plot = FALSE)
      
      conc <- (1000/0.91 * dpcr$lambda)*(vol.mix/vol.sample)
      conc_lower <- (1000/0.91 * dpcr$lower)*(vol.mix/vol.sample)
      conc_upper <- (1000/0.91 * dpcr$upper)*(vol.mix/vol.sample)
      
    } else{ 
      conc_lower <- 0
      conc_upper <- 0
      conc <- 0
    }
    

  }
}
