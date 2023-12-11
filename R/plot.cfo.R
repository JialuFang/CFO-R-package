#' Plot the simulation results for CFO-type and aCFO-type designs
#'
#' Plot the objects returned by other functions, including (1) dose allocation of a single trial;
#' (2) operating characteristics of multiple simulations, including selesction percentage and the
#' number of patients treated at each dose
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
#' @author Jialu Fang
#' 
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics axis barplot mtext par plot
#' @import ggplot2
#' @export
#'
#' @examples
#' ############design without late-onset outcomes################
#' nsimu <- 100; phi <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' add.args=list(alp.prior=phi, bet.prior=1-phi)
#' ## CFO design
#' CFOtrial <- CFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args, accumulation = FALSE)
#' plot.cfo(CFOtrial)
#' CFOsimu <- CFO.oc (nsimu, design='CFO', phi, p.true, ncohort, init.level, cohortsize,
#'                     tau=NaN, accrual=NaN, tite.dist=NaN, accrual.dist=NaN, add.args)
#' plot.cfo(CFOsimu)
#' ## aCFO design
#' aCFOtrial <- CFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args,
#'                     accumulation = TRUE)
#' plot.cfo(aCFOtrial)
#' aCFOsimu <- CFO.oc (nsimu, design='aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                     tau=NaN, accrual=NaN, tite.dist=NaN, accrual.dist=NaN, add.args)
#' plot.cfo(aCFOsimu)
#' 
#' 
#' ##############design with late-onset outcomes################
#' tau <- 3; accrual <- 6; tite.dist <- 2; accrual.dist <- 1
#' ## TITE-CFO design
#' TITECFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist,
#'                 accrual.dist, design='TITE-CFO', init.level, add.args)
#' plot.cfo(TITECFOtrial)
#' TITECFOsimu <- CFO.oc (nsimu, design='TITE-CFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot.cfo(TITECFOsimu)
#' ## TITE-aCFO design
#' TITEaCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist,
#'                  accrual.dist,design='TITE-aCFO', init.level, add.args)
#' plot.cfo(TITEaCFOtrial)
#' TITEaCFOsimu <- CFO.oc (nsimu, design='TITE-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot.cfo(TITEaCFOsimu)
#' ## fCFO design
#' fCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist,
#'                 accrual.dist, design='fCFO', init.level, add.args)
#' plot.cfo(fCFOtrial)
#' fCFOsimu <- CFO.oc (nsimu, design='fCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot.cfo(fCFOsimu)
#' ## f-aCFO design
#' faCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist,
#'                 accrual.dist, design='f-aCFO', init.level, add.args)
#' plot.cfo(faCFOtrial)
#' faCFOsimu <- CFO.oc (nsimu, design='f-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot.cfo(faCFOsimu)
#' 
#' ##############two-dim design with late-onset outcomes################
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#' 
#' ## plot the single trail returned by CFO2d.sim()
#' CFO2dtrail <- CFO2d.sim(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3)
#' plot.cfo(CFO2dtrail)
#' 
#' ## plot the multiple trails returned by CFO2d.oc()
#' CFO2dsim <- CFO2d.oc(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3, n_sim = 100)
#' plot.cfo(CFO2dsim)
#' 
plot.cfo<- function (x,..., name = deparse(substitute(x)))
{
  new.obj = unlist(strsplit(name, split = "\\$"))
  strpattern = "none"
  if (length(new.obj) >= 2) {
    strpattern = new.obj[2]
  }
  assign("objectPlot", get(new.obj[1]))
  if (!is.element(strpattern, c("none", names(objectPlot)))) {
    warning("Please double check and specify the variable to be plotted...\n")
  }
  else {
    if (!is.null(objectPlot$simu.oc)) { #plot for one-dim multiple simulations
      message("1d oc")
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      dev.flush()
      dev.hold()
      par(mar = c(5, 6, 4, 2))
      bplot = barplot(objectPlot$selPercent*100, ylab = "selection percentage (%)",
                      ylim = c(0, 100), cex.names = 1, xaxt = "n",
                      cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$selPercent)))
      mtext("Selection percentage", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
      dev.flush()
      dev.hold()
      bplot = barplot(objectPlot$nPatients, ylab = "number of patients",
                      ylim = c(0, max(objectPlot$nPatients))*1.3, cex.names = 1,
                      beside = FALSE, xaxt = "n", cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$nPatients)))
      mtext("Patient allocation", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
      dev.flush()
      dev.hold()
      bplot = barplot(objectPlot$nTox, ylab = "number of toxicities",
                      ylim = c(0, max(objectPlot$nTox)*1.3), cex.names = 1,
                      beside = FALSE, xaxt = "n", cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$nTox)))
      mtext("Observed toxicity", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
    } else if (!is.null(objectPlot$dose.list)) { #plot for one-dim single trial
      if (!is.null(objectPlot$total.time)){ #plot for one-dim single trial (design with lateonset outcomes)
        dose <- objectPlot$dose.list
        DLT <- objectPlot$DLT.list
        ncohort <- length(objectPlot$dose.list)
        cohortsize <- sum(objectPlot$dose.ns)/ncohort
        
        # Generate y_labels
        y_labels <- seq(1, max(dose))
        
        # Generate sequences for each patient
        sequences <- objectPlot$enter.times
        
        # Generate dose_levels for each patient
        dose_levels <- rep(dose, each = cohortsize)
        
        # Generate DLT_observed for each patient
        DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
        
        new_seq <- ifelse(objectPlot$DLT.times!=0, sequences+objectPlot$DLT.times, NaN)
        new_y <- ifelse(objectPlot$DLT.times!=0, dose_levels+0.3, NaN)
        
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
      else{  #plot for one-dim single trial (design without lateonset outcomes)
        dose <- objectPlot$dose.list
        DLT <- objectPlot$DLT.list
        ncohort <- length(objectPlot$dose.list)
        cohortsize <- sum(objectPlot$dose.ns)/ncohort
        
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
    } else if (!is.null(objectPlot$p.true)){ #plot for two-dim simulations
      if (is.null(objectPlot$selPercent)){
        dose <- objectPlot$dose
        DLT <- objectPlot$DLT
        ncohort <- length(objectPlot$DLT)
        cohortsize <- sum(objectPlot$dose.ns)/ncohort
        dim <- dim(objectPlot$p.true)
        
        # Generate y_labels
        y_labels <- expand.grid(1:dim[1], 1:dim[2])
        y_labels <- apply(y_labels, 1, function(x) paste('(', x[1], ',', x[2], ')'))
        
        # Generate sequences for each patient
        sequences <- 1:(ncohort * cohortsize)
        
        # Generate dose_levels for each patient
        dose_levels <- rep(match(apply(dose, 1, function(x) paste('(', x[1], ',', x[2], ')')), y_labels), each = cohortsize)
        
        # Generate DLT_observed for each patient
        DLT_observed <- unlist(mapply(function(dlt, size) c(rep(1, dlt), rep(0, size - dlt)),DLT, rep(cohortsize, ncohort)))
        
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
      } else {
        attributesToPlot <- c("selPercent", "dose.ns", "DLT.ns")
        titles <- c("MTD selection rate", "Average patients allocation", "Average DLT observed")
        ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs")
        
        par(mfrow = c(3, 1))
        
        # Loop through each attribute and create a plot
        for (i in seq_along(attributesToPlot)) {
          attr <- attributesToPlot[i]
          # Check if the attribute exists in the objectPlot
          if (!is.null(objectPlot[[attr]])) {
            # Extract the matrix
            matrixToPlot <- objectPlot[[attr]]
            
            # Convert the matrix to a vector by column
            matrixVector <- as.vector(matrixToPlot)
            
            # Convert to percentages only for selPercent
            if (attr == "selPercent") {
              matrixVector <- matrixVector * 100
            }
            
            # Create x-axis labels
            dimMatrix <- dim(matrixToPlot)
            xLabels <- expand.grid(row = 1:dimMatrix[1], col = 1:dimMatrix[2])
            xLabels <- apply(xLabels, 1, function(x) paste("(", x[1], ",", x[2], ")", sep = ""))
            
            # Create the bar plot with horizontal x-axis labels
            barplot(matrixVector, names.arg = xLabels, las = 2,
                    xlab = "Index (row,col)", ylab = ylabels[i],
                    main = titles[i])
          }
        }
        par(mfrow = c(1, 1))
      }
    } 
  }
}



