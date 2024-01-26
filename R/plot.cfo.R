#' Plot the simulation results for CFO-type and aCFO-type designs
#'
#' Plot the objects returned by other functions, including (1) dose allocation of a single trial;
#' (2) operating characteristics of multiple simulations, including MTD selesction percentage, 
#' the averaged number of patients allocated for different doses in one simulation and the averaged 
#' number of DLT observed for different doses in one simulation.
#' 
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User doesn't need to input this parameter.
#'
#' @details \code{plot()} returns a figure or a series of figures depending on the object entered.
#'           Additionally, in the example, we set \code{nsimu=100} for testing time considerations. 
#'           In reality, \code{nsimu} is typically set to 5000 to ensure the accuracy of the results.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#' 
#' @author Jialu Fang and Wenliang Wang
#' 
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics axis barplot mtext par plot
#' @importFrom graphics abline arrows legend points
#' @import ggplot2
#' @export
#'
#' @examples
#' ############design without late-onset outcomes################
#' nsimu <- 100; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <- c(0.01, 0.02, 0.05, 0.20, 0.30, 0.50, 0.70)
#' target <- 0.2; prior.para = list(alp.prior = target, bet.prior = 1 - target)
#' ## CFO design
#' CFOtrial <- CFO.simu(target, p.true, ncohort, init.level = 1, cohortsize = 3, design = 'CFO',
#'             prior.para)
#' plot(CFOtrial)
#' CFOoc <- CFO.oc (nsimu, design='CFO', target, p.true, ncohort, init.level, cohortsize,
#'                     tau = NaN, tte.para = NaN, accrual = NaN, accrual.dist = NaN, prior.para)
#' plot(CFOoc)
#' ## aCFO design
#' aCFOtrial <- CFO.simu(target, p.true, ncohort, init.level = 1, cohortsize = 3, design = 'aCFO', 
#'               prior.para)
#' plot(aCFOtrial)
#' aCFOoc <- CFO.oc (nsimu, design = 'aCFO', target, p.true, ncohort, init.level, cohortsize,
#'                     tau = NaN, tte.para = NaN, accrual = NaN, accrual.dist = NaN, prior.para)
#' plot(aCFOoc)
#' 
#' 
#' ##############design with late-onset outcomes################
#' tau <- 3; accrual <- 2; tte.para <- 0.5; accrual.dist <- 'unif'
#' ## TITE-CFO design
#'  TITECFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                  accrual.dist, design = 'TITE-CFO', init.level = 1, 
#'                  prior.para = list(alp.prior = target, bet.prior = 1 - target))
#' plot(TITECFOtrial)
#' TITECFOoc <- CFO.oc (nsimu, design='TITE-CFO', target, p.true, ncohort, init.level, cohortsize,
#'                       tau, tte.para, accrual, accrual.dist, prior.para)
#' plot(TITECFOoc)
#' ## TITE-aCFO design
#' TITEaCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'TITE-aCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target))
#' plot(TITEaCFOtrial)
#' TITEaCFOoc <- CFO.oc (nsimu, design='TITE-aCFO', target, p.true, ncohort, init.level, cohortsize,
#'                       tau, tte.para, accrual, accrual.dist, prior.para)
#' plot(TITEaCFOoc)
#' ## fCFO design
#' fCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'fCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target))
#' plot(fCFOtrial)
#' fCFOoc <- CFO.oc (nsimu, design = 'fCFO', target, p.true, ncohort, init.level, cohortsize,
#'                       tau, tte.para, accrual, accrual.dist, prior.para)
#' plot(fCFOoc)
#' ## f-aCFO design
#' faCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'f-aCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target))
#' plot(faCFOtrial)
#' faCFOoc <- CFO.oc (nsimu, design = 'f-aCFO', target, p.true, ncohort, init.level, cohortsize,
#'                       tau, tte.para, accrual, accrual.dist, prior.para)
#' plot(faCFOoc)
#' 
#' ##############two-dim design with late-onset outcomes################
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#' 
#' ## plot the single trail returned by CFO2d.simu()
#' CFO2dtrail <- CFO2d.simu(target = 0.3, p.true = p.true, ncohort = 20, cohortsize = 3)
#' plot(CFO2dtrail)
#' 
#' ## plot the multiple trails returned by CFO2d.oc()
#' CFO2doc <- CFO2d.oc(target = 0.3, p.true = p.true, ncohort = 20, cohortsize = 3, 
#'             seeds = 1:100, nsimu = 100)
#' plot(CFO2doc)
 
plot.cfo<- function (x,..., name = deparse(substitute(x)))
{
  new.obj = unlist(strsplit(name, split = "\\$"))
  strpattern = "none"
  if (length(new.obj) >= 2){
    strpattern = new.obj[2]
  }
  assign("objectPlot", get(new.obj[1]))
  if (!is.element(strpattern, c("none", names(objectPlot)))) {
    warning("Please double check and specify the variable to be plotted...\n")
  }
  else {
    ###############################################################################
    ############################plot for CFO.oc()###############################
    ###############################################################################
    if (!is.null(objectPlot$simu.setup)) { #plot for one-dim multiple simulations
      if(is.null(dim(objectPlot$selpercent))){
        attributesToPlot <- c("selpercent", "npatients", "ntox")
        titles <- c("MTD selection rate", "Average patients allocation", "Average DLT observed")
        ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs")
        
        par(mfrow = c(3, 1))
        
        # Loop through each attribute and create a plot
        for (i in seq_along(attributesToPlot)) {
          attr <- attributesToPlot[i]
          # Check if the attribute exists in the objectPlot
          if (!is.null(objectPlot[[attr]])) {
            # Extract the vector
            vectorToPlot <- objectPlot[[attr]]
            
            # Convert to percentages only for selPercent
            if (attr == "selpercent") {
              vectorToPlot <- vectorToPlot * 100
            }
            
            # Create the bar plot with horizontal x-axis labels
            bplot <- barplot(vectorToPlot, ylab = ylabels[i], main = titles[i], xlab = "Dose level",
                             cex.names = 1, xaxt = "n", ylim = c(0, max(vectorToPlot))*1.3,
                             cex.lab = 1.3)
            axis(1, at = bplot, labels = seq(1, length(objectPlot[[attr]])))
          }
        }
        par(mfrow = c(1, 1))
      }
      else if(length(dim(objectPlot$selPercent))==2) {
        #   attributesToPlot <- c("selPercent", "dose.ns", "DLT.ns")
        #   titles <- c("MTD selection rate", "Average patients allocation", "Average DLT observed")
        #   ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs")
        #   
        #   par(mfrow = c(3, 1))
        #   
        #   # Loop through each attribute and create a plot
        #   for (i in seq_along(attributesToPlot)) {
        #     attr <- attributesToPlot[i]
        #     # Check if the attribute exists in the objectPlot
        #     if (!is.null(objectPlot[[attr]])) {
        #       # Extract the matrix
        #       matrixToPlot <- objectPlot[[attr]]
        #       
        #       # Convert the matrix to a vector by column
        #       matrixVector <- as.vector(matrixToPlot)
        #       
        #       # Convert to percentages only for selPercent
        #       if (attr == "selPercent") {
        #         matrixVector <- matrixVector * 100
        #       }
        #       
        #       # Create x-axis labels
        #       dimMatrix <- dim(matrixToPlot)
        #       xLabels <- expand.grid(row = 1:dimMatrix[1], col = 1:dimMatrix[2])
        #       xLabels <- apply(xLabels, 1, function(x) paste("(", x[1], ",", x[2], ")", sep = ""))
        #       
        #       # Create the bar plot with horizontal x-axis labels
        #       barplot(matrixVector, names.arg = xLabels, las = 2,
        #               xlab = "Combined dose level", ylab = ylabels[i], main = titles[i])
        #     }
        #   }
        #   par(mfrow = c(1, 1))
      }
    }
    
    
    ###############################################################################
    #########################plot for XXX.simu()###################################
    ###############################################################################
    if (!is.null(objectPlot$correct)) { 
      if(length(objectPlot$MTD) == 1){
        if (!is.null(objectPlot$totaltime)){ #plot for CFO.simu()
          dose <- objectPlot$cohortdose
          DLT <- objectPlot$patientDLT
          ncohort <- length(objectPlot$cohortdose)
          cohortsize <- sum(objectPlot$npatients)/ncohort
          
          # Generate y_labels
          y_labels <- seq(1, max(dose))
          
          # Generate sequences for each patient
          sequences <- objectPlot$entertimes
          
          # Generate dose_levels for each patient
          dose_levels <- rep(dose, each = cohortsize)
          
          # Generate DLT_observed for each patient
          DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
          
          new_seq <- ifelse(objectPlot$DLTtimes!=0, sequences+objectPlot$DLTtimes, NaN)
          new_y <- ifelse(objectPlot$DLTtimes!=0, dose_levels+0.3, NaN)
          
          df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
          dfnew <- data.frame(sequence = sequences, dose_levels = dose_levels, new_seq = new_seq, new_y = new_y)
          dfnew <- na.omit(dfnew)
          
          # Create the plot
          p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
            geom_point(aes(shape = factor(DLT_observed,levels=c(0,1,2))), color = 'black', size = 2.5) +
            geom_step(direction = 'hv', color = 'black') +
            scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
            labs(x = "Time (in months)", 
                 y = "Dose level",
                 fill = 'DLT observed') +
            theme_minimal() +
            theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
            scale_shape_manual(values = c(1, 16, 4), labels = c('DLT not observed', 'DLT observed',"DLT time"), drop = FALSE)
          
          for (row in 1:(nrow(dfnew))){
            xuse=c(dfnew[row,"sequence"],dfnew[row,"new_seq"])
            yuse=c(dfnew[row,"dose_levels"],dfnew[row,"new_y"])
            dfuse <-data.frame(xuse=xuse, yuse=yuse)
            p <- p + geom_point(aes(x = xuse[2], y = yuse[2]), shape = 4,size = 2.5, data = dfuse)+
              geom_step(aes(x = xuse, y = yuse), data = dfuse,direction = 'vh',
                        linetype = 2)
          }
          print(p)
        }
        else{  #plot for lateonset.simu()
          dose <- objectPlot$cohortdose
          DLT <- objectPlot$patientDLT
          ncohort <- length(objectPlot$cohortdose)
          cohortsize <- sum(objectPlot$npatients)/ncohort
          
          # Generate y_labels
          y_labels <- seq(1, max(dose))
          
          # Generate sequences for each patient
          sequences <- 1:(ncohort * cohortsize)
          
          # Generate dose_levels for each patient
          dose_levels <- rep(dose, each = cohortsize)
          
          # Generate DLT_observed for each patient
          DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
          
          df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
          
          # Create the plot
          p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
            geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
            geom_step(direction = 'hv', color = 'black') +
            scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
            labs(x = "Sequence of patients treated", 
                 y = "Dose level",
                 fill = 'DLT observed') +
            theme_minimal() +
            theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
            scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
          
          # Display the plot
          print(p)
        }
      }
      else{
        dose <- objectPlot$cohortdose
        DLT <- objectPlot$patientDLT   ###need to change!!!!
        ncohort <- dim(objectPlot$cohortdose)[1]
        cohortsize <- sum(objectPlot$npatients)/ncohort
        dim <- dim(objectPlot$ntox)
        
        # Generate y_labels
        y_labels <- expand.grid(1:dim[1], 1:dim[2])
        y_labels <- apply(y_labels, 1, function(x) paste('(', x[1], ',', x[2], ')'))
        
        # Generate sequences for each patient
        sequences <- 1:(ncohort * cohortsize)
        
        # Generate dose_levels for each patient
        dose_levels <- rep(match(apply(dose, 1, function(x) paste('(', x[1], ',', x[2], ')')), y_labels), each = cohortsize)
        
        # Generate DLT_observed for each patient
        DLT_observed <- unlist(mapply(function(dlt, size) c(rep(1, dlt), rep(0, size - dlt)), DLT, rep(cohortsize, ncohort)))
        
        df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
        
        # Create the plot
        p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
          geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
          geom_step(direction = 'hv', color = 'black') +
          scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
          labs(x = "Sequence of patients treated", 
               y = "Combined dose level",
               fill = 'DLT observed') +
          theme_minimal() +
          theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
          scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
        # Display the plot
        print(p)
      }
    }
    
    ###############################################################################
    #########################plot for select.mtd()###################################
    ###############################################################################
    if (!is.null(objectPlot$p_est)){
      
      if (objectPlot$MTD[1] == 99) {
        warning("All tested doses are overly toxic. No MTD is selected!\n")
      }
      else {
        if (!is.null(objectPlot$p_est)) {
          
          par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
          if (length(objectPlot$MTD) >= 2) {
            p_est.comb=objectPlot$p_est
            rownames(p_est.comb)=1:dim(p_est.comb)[1]
            colnames(p_est.comb)=1:dim(p_est.comb)[2]
            barplot(p_est.comb,beside=TRUE,ylab="DLT rate",
                    ylim=c(0,round(max(p_est.comb,na.rm=TRUE)*1.5,1)),xlab="Drug B",legend.text=rownames(p_est.comb),
                    args.legend=list(title="Drug A",horiz=TRUE,x="top"))
          }
          else {
            p_est = objectPlot$p_est
            p_hat = p_est[, 2]
            ci = p_est[, 3]
            ci = gsub("[\\(\\)]", "", ci)
            conf.intv = matrix(unlist(strsplit(ci, ",")),
                               byrow = TRUE, ncol = 2)
            if (p_est[1, 2] == "----") {
              warning("The trial is stopped since the lowest dose is too toxic.\n")
            }
            else {
              numbs = ifelse(sum(p_hat == "----") ==
                               0, length(p_hat), min(which(p_hat ==
                                                             "----")) - 1)
              numbs2 = length(p_hat)
              phatx = as.numeric(as.character(p_hat[1:numbs]))
              lwr = as.numeric(as.character(conf.intv[1:numbs,
                                                      1]))
              upr = as.numeric(as.character(conf.intv[1:numbs,
                                                      2]))
              par(mar = c(5, 5, 4, 2))
              plot(1:numbs2, ylim = c(0, 1), xlab = "Dose level",
                   ylab = "DLT rate", pch = "", xaxt = "n",
                   cex.lab = 1.3)
              axis(1, at = 1:numbs2, labels = 1:numbs2)
              abline(h = objectPlot$target, lty = 2,
                     col = 2)
              points(1:numbs, phatx, pch = 19)
              arrows(x0 = 1:numbs, x1 = 1:numbs, y0 = lwr,
                     y1 = upr, code = 3, angle = 90, length = 0.1)
              if (numbs < numbs2) {
                points((numbs + 1):numbs2, seq(min(1,
                                                   max(phatx, na.rm = T) + 0.05), min(max(phatx,
                                                                                          na.rm = T) + 0.2, 1), length = numbs2 -
                                                 numbs), pch = "*", cex = 1.5)
                legend("topleft", "*   no patient treated")
              }
            }
          }
        }
      }

  
    }
  }
}

