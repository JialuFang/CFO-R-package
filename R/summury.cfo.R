#' Generate descriptive summary for objects returned by other functions
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param object the object returned by other functions.
#' @param ... ignored arguments
#'
#'
#' @details \code{summary()} prints the objects returned by other functions. Additionally, in the example, 
#'          we set \code{nsimu=100} for testing time considerations. In reality, \code{nsimu} is typically 
#'          set to 5000 to ensure the accuracy of the results.
#'
#' @return \code{summary()} prints the objects returned by other functions.
#'
#' @author Jialu Fang and Wenliang Wang
#'
#' @export
#' 
#' @examples
#' ## settings for 1dCFO
#' nsimu <- 100; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <- c(0.02, 0.05, 0.20, 0.28, 0.34, 0.40, 0.44)
#' target <- 0.2; prior.para=list(alp.prior=target, bet.prior=1-target)
#' tau <- 3; acrrual.rate <- 2; tte.para <- 0.5; accrual.dist <- 'unif'
#' 
#' ## summarize the object returned by CFO.next()
#' decision <- CFO.next(target = 0.2, cys = c(0, 1, 0), cns = c(3, 6, 0), currdose = 3, prior.para)
#' summary(decision)
#' 
#' ## summarize the object returned by CFO.simu()
#' 
#' aCFOtrial <- CFO.simu(target, p.true, ncohort, init.level = 1, cohortsize = 3, design = 'aCFO',
#'                       prior.para = list(alp.prior = target, bet.prior = 1 - target))
#' summary(aCFOtrial)
#' 
#' ## summarize the object returned by lateonset.next()
#' enter.times<- c(0, 0.266, 0.638, 1.54, 2.48, 3.14, 3.32, 4.01, 4.39, 5.38, 5.76,
#'                6.54, 6.66, 6.93, 7.32, 7.65, 8.14, 8.74)
#' dlt.times<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0.995, 0, 0, 0, 0, 0, 0, 0, 2.58)
#' current.t<- 9.41
#' doses<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4)
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'f-aCFO', tau = 3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' 
#' ## summarize the object returned by lateonset.simu()
#' faCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, acrrual.rate,
#'                 accrual.dist, design = 'f-aCFO', init.level, prior.para)
#' summary(faCFOtrial)
#' 
#' ## summarize the object returned by CFO.oc()
#' faCFOoc <- CFO.oc (nsimu, design = 'f-aCFO', target, p.true, ncohort, init.level, cohortsize,
#'                       tau, tte.para, acrrual.rate, accrual.dist, prior.para)
#' summary(faCFOoc)
#' 
#' ## settings for 2dCFO
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#' 0.10, 0.15, 0.30, 0.45, 0.55,
#' 0.15, 0.30, 0.45, 0.50, 0.60), 
#' nrow = 3, ncol = 5, byrow = TRUE)
#' 
#' cns <- matrix(c(3, 3, 0,
#'                 0, 6, 0,
#'                 0, 0, 0), 
#'               nrow = 3, ncol = 3, byrow = TRUE)
#' cys <- matrix(c(0, 1, 0,
#'                 0, 2, 0,
#'                 0, 0, 0), 
#'               nrow = 3, ncol = 3, byrow = TRUE)
#' currdose <- c(2,3)
#' 
#' ## summarize the object returned by CFO2d.next()
#' decision <- CFO2d.next(0.3, cys, cns, currdose = currdose)
#' summary(decision)
#' 
#' ## summarize the object returned by CFO2d.simu()
#' CFO2dtrail <- CFO2d.simu(target=0.3, p.true=p.true, ncohort = 20, cohortsize = 3)
#' summary(CFO2dtrail)
#' 
#' ## summarize the object returned by CFO2d.oc()
#' CFO2doc <- CFO2d.oc(target=0.3, p.true=p.true, ncohort = 20, cohortsize = 3, nsimu = 100)
#' summary(CFO2doc)
#' 
summary.cfo<- function (object, ...)
{
  
  ###############################################################################
  ############################summary for XXX.oc()###############################
  ###############################################################################
  
  if (!is.null(object$simu.setup)) {
    if(is.null(dim(object$selpercent))){
      if (object$percentstop == 0){
        cat("No instance of early stopping was observed in",
            object$simu.setup$nsimu, "simulations. \n")
      }else{
        nstop = object$percentstop*object$simu.setup$nsimu
        cat("In", object$simu.setup$nsimu, "simulations, early stopping occurred",
            nstop, "times \n")
        cat("Among simulations where early stopping did not occur: \n")
      }
      
      cat("Selection percentage at each dose level:\n")
      cat(formatC(object$selpercent, digits = 3, format = "f"),
          sep = "  ", "\n")
      cat("Average number of patients treated at each dose level:\n")
      cat(formatC(object$npatients, digits = 3, format = "f"),
          sep = "  ", "\n")
      cat("Average number of toxicity observed at each dose level:\n")
      cat(formatC(object$ntox,  digits = 3, format = "f"),
          sep = "  ", "\n")
      cat("Percentage of correct selection of the MTD:", 
          formatC(object$MTDsel, digits = 3, format = "f"), "\n")
      cat("Percentage of patients allocated to the MTD:", 
          formatC(object$MTDallo, digits = 3, format = "f"), "\n")
      cat("Percentage of selecting a dose above the MTD:",
          formatC(object$oversel, digits = 3, format = "f")," \n")
      cat("Percentage of allocating patients at dose levels above the MTD:",
          formatC(object$overallo, digits = 3, format = "f")," \n")
      cat("Percentage of the patients suffering DLT:",
          formatC(object$averDLT, digits = 3, format = "f")," \n")
      
      if (!is.null(object$averdur)){
        cat("Average trial duration:",
            formatC(object$averdur, digits = 3, format = "f")," \n")
      }
    }
    else if(length(dim(object$selPercent))==2) {
      # Summary for 2dCFO multiple trail simulation
      cat("Selection percentage at each dose level:\n")
      print(object$selPercent)
      cat("Average number of patients treated at each dose level:\n")
      print(object$dose.ns)
      cat("Average number of toxicity observed at each dose level:\n")
      print(object$DLT.ns)
      cat("Percentage of correct selection of the MTD:", 
          formatC(object$pcorrect, digits = 3, format = "f"), "\n")
      cat("Percentage of patients allocated to the MTD:", 
          formatC(object$npercent, digits = 3, format = "f"), "\n")
      cat("Percentage of allocating patients at dose levels above the MTD:",
          formatC(object$ptoxic, digits = 3, format = "f")," \n")
      cat("Percentage of the patients suffering DLT:",
          formatC(object$ntox/sum(object$dose.ns), digits = 3, format = "f")," \n")
    }
  }
  
  ###############################################################################
  #########################summary for XXX.simu()###############################
  ###############################################################################
  
  if(!is.null(object$correct)){ ###summary for XXX.simu()
    if (length(object$MTD) == 1) {  ###summary for one-dim XXX.simu()
      if (object$MTD == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        cat("The selected MTD is dose level", object$MTD, "\n")
        cat("For",length(object$cohortdose),"cohorts, the dose level assigned to each cohort is: \n")
        cat(formatC(object$cohortdose, format = "d"), sep = "  ", "\n")
        cat("Number of toxicity observed at each dose level:\n")
        cat(formatC(object$ntox, format = "d"), sep = "  ", "\n")
        cat("Number of patients treated at each dose level:\n")
        cat(formatC(object$npatients, format = "d"), sep = "  ", "\n")
        if (!is.null(object$totaltime)){
          cat("The duration of the trial in months:",
              formatC(object$totaltime, digits = 3, format = "f")," \n")
        }
      }
    } else {  ###summary for two-dim XXX.simu()
      if (object$MTD[1] == 99 | object$MTD[2] == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        # Summary for 2dCFO single trail simulation
        cat("The selected MTD is dose level (", object$MTD[1], ",",object$MTD[2], ")\n\n")
      }
      # print assgined dosage for each cohort
      doses <- object$cohortdose
      cohort_data <- data.frame(
        cohort = 1:nrow(doses),
        dose_A = doses[, 1],
        dose_B = doses[, 2]
      )
      print(cohort_data, row.names = FALSE)
      cat("\n")
      cat("Number of toxicity observed at each dose level:\n")
      print(object$ntox)
      cat("\n")
      cat("Number of patients treated at each dose level:\n")
      print(object$npatients)
    }
  }
  
  ###############################################################################
  #########################summary for XXX.next()################################
  ###############################################################################
  
  if(!is.null(object$decision)){
    if(length(object$decision)==2){ ##summary for two dim XXX.next()
      if (is.na(object$overtox)) {
        cat("All tested doses are not overly toxic. \n\n")
      } else {
        cat("Dose level", object$overtox, "and all levels above exhibit excessive toxicity", "\n")
      }
      cat("The decision regarding the direction of movement for drug A is", object$decision[1], "\n")
      cat("The decision regarding the direction of movement for drug B is", object$decision[2], "\n")
      cat("The next cohort will be assigned to dose level (", object$nextdose[1],",",object$nextdose[2],")", "\n")
    } else {                       ##summary for one dim XXX.next()
      if (is.na(object$overtox)) {
        cat("All tested doses are not overly toxic \n\n")
      } else {
        cat("Dose level", object$overtox, "and all levels above exhibit excessive toxicity", "\n")
      }
      if (object$decision == "stop"){
        cat("The lowest dose level is overly toxic. We terminate the entire trial for safety.")
      }else{
        cat("The current dose level is", object$currdose, "\n")
        cat("The decision regarding the direction of movement is", object$decision, "\n")
        cat("The next cohort will be assigned to dose level", object$nextdose, "\n")
      }
    }
  }
  
  ###############################################################################
  #########################summary for CFO.selectmtd()##############################
  ###############################################################################
  
  if (!is.null(object$p_est)){ ##summary for CFO.selectmtd()
    if (length(object$MTD) == 1) { ##summary for one dim CFO.selectmtd()
      if (object$MTD == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        cat("The MTD is dose level ", object$MTD, "\n\n")
      }
      cat("Dose    Posterior DLT             95%                  \n",
          sep = "")
      cat("Level     Estimate         Credible Interval   Pr(toxicity>",
          object$target, "|data)\n", sep = "")
      for (i in 1:nrow(object$p_est)) {
        cat(" ", i, "        ", as.character(object$p_est[i,2]), "         ", 
            as.character(object$p_est[i,3]), "         ", as.character(object$p_overdose[i]),
            "\n")
      }
      cat("NOTE: no estimate is provided for the doses at which no patient was treated.\n")
    }
  }
}



