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
#' @examples
#' ## settings for 1dCFO
#' nsimu <- 100; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' phi <- 0.2; add.args=list(alp.prior=phi, bet.prior=1-phi)
#' tau <- 3; accrual <- 6; tite.dist <- 2; accrual.dist <- 1
#' 
#' ## summarize the object returned by CFO.next()
#' decision <- CFO.next(phi=0.2, cys=c(0,1,0), cns=c(3,6,0), curDose=3, add.args)
#' summary(decision)
#' 
#' ## summarize the object returned by CFO.simu()
#' 
#' aCFOtrial <- CFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3, design='aCFO',
#'                       add.args=list(alp.prior=phi, bet.prior=1-phi))
#' summary(aCFOtrial)
#' 
#' ## summarize the object returned by lateonset.next()
#' enter.times<-c(0,0.082,0.343,0.554,1.18,1.88,2.15,2.68,2.74,3.30,3.99,
#'                4.52,5.51,5.93,6.18,6.68,7.05,7.92)
#' dlt.times<-c(0,0,0,0,0,0,0,0,0,0.934,0,0,0,0,0,1.5,0.6962,0)
#' current.t<-8.134413
#' doses<-c(1,1,1,2,2,2,3,3,3,4,4,4,3,3,3,4,4,4)
#' decision <- lateonset.next(phi, p.true, curDose=4, design='fCFO', tau=3, enter.times, dlt.times,  
#'                            current.t, doses, add.args)
#' summary(decision)
#' 
#' ## summarize the object returned by lateonset.simu()
#' faCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, 
#'                 accrual.dist, design='f-aCFO', init.level, add.args)
#' summary(faCFOtrial)
#' 
#' ## summarize the object returned by CFO.oc()
#' faCFOsimu <- CFO.oc (nsimu, design='f-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' summary(faCFOsimu)
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
#' curDose <- c(2,3)
#' 
#' ## summarize the object returned by CFO2d.next()
#' decision <- CFO2d.next(0.3, cys, cns, curDose = curDose)
#' summary(decision)
#' 
#' ## summarize the object returned by CFO2d.sim()
#' CFO2dtrail <- CFO2d.simu(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3)
#' summary(CFO2dtrail)
#' 
#' ## summarize the object returned by CFO2d.oc()
#' CFO2dsim <- CFO2d.oc(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3, nsimu = 100)
#' summary(CFO2dsim)
#' 
summary.cfo<- function (object, ...)
{
  if (!is.null(object$simu.oc)) {
    if (object$simu.oc$errStop == 0){
      cat("No instance of early stopping was observed in",
          object$simu.setup$nsimu, "simulations. \n")
    }else{
      cat("In", object$simu.setup$nsimu, "simulations, early stopping occurred",
          object$simu.oc$errStop, "times \n")
      cat("Among simulations where early stopping did not occur: \n")
    }
    
    cat("Selection percentage at each dose level:\n")
    cat(formatC(object$selPercent, digits = 3, format = "f"),
        sep = "  ", "\n")
    cat("Average number of patients treated at each dose level:\n")
    cat(formatC(object$dose.ns, digits = 3, format = "f"),
        sep = "  ", "\n")
    cat("Average number of toxicity observed at each dose level:\n")
    cat(formatC(object$DLT.ns,  digits = 3, format = "f"),
        sep = "  ", "\n")
    cat("Percentage of correct selection of the MTD:", 
        formatC(object$simu.oc$MTDSel, digits = 3, format = "f"), "\n")
    cat("Percentage of patients allocated to the MTD:", 
        formatC(object$simu.oc$MTDAllo, digits = 3, format = "f"), "\n")
    cat("Percentage of selecting a dose above the MTD:",
        formatC(object$simu.oc$overSel, digits = 3, format = "f")," \n")
    cat("Percentage of allocating patients at dose levels above the MTD:",
        formatC(object$simu.oc$overAllo, digits = 3, format = "f")," \n")
    cat("Percentage of the patients suffering DLT:",
        formatC(object$simu.oc$averDLT, digits = 3, format = "f")," \n")
    
    if (!is.null(object$simu.oc$averDur)){
      cat("Average trial duration:",
          formatC(object$simu.oc$averDur, digits = 3, format = "f")," \n")
    }
  }
  
  if(length(dim(object$selPercent))==2) {
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
  
  if(!is.null(object$MTD)){
    if (length(object$MTD) == 1) {
      if (object$MTD == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        cat("The selected MTD is dose level", object$MTD, "\n")
        cat("For",length(object$dose.list),"cohorts, the dose level assigned to each cohort is: \n")
        cat(formatC(object$dose.list, format = "d"), sep = "  ", "\n")
        cat("Number of toxicity observed at each dose level:\n")
        cat(formatC(object$DLT.ns, format = "d"), sep = "  ", "\n")
        cat("Number of patients treated at each dose level:\n")
        cat(formatC(object$dose.ns, format = "d"), sep = "  ", "\n")
        if (!is.null(object$total.time)){
          cat("The duration of the trial in months:",
              formatC(object$total.time, digits = 3, format = "f")," \n")
        }
      }
    } else {
      if (object$MTD[1] == 99 | object$MTD[2] == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        # Summary for 2dCFO single trail simulation
        cat("The selected MTD is dose level (", object$MTD[1], ",",object$MTD[2], ")\n\n")
      }
      # print assgined dosage for each cohort
      doses <- object$dose
      cohort_data <- data.frame(
        cohort = 1:nrow(doses),
        dose_A = doses[, 1],
        dose_B = doses[, 2]
      )
      print(cohort_data, row.names = FALSE)
      cat("\n")
      cat("Number of toxicity observed at each dose level:\n")
      print(object$DLT.ns)
      cat("\n")
      cat("Number of patients treated at each dose level:\n")
      print(object$dose.ns)
    }
  }
  
  if(!is.null(object$decision)){
    if(length(object$decision)==2){
      if (is.na(object$overTox)) {
        cat("All tested doses are not overly toxic. \n\n")
      } else {
        cat("Dose level", object$overTox, "and all levels above exhibit excessive toxicity", "\n")
      }
      cat("The decision regarding the direction of movement for drug A is", object$decision[1], "\n")
      cat("The decision regarding the direction of movement for drug B is", object$decision[2], "\n")
      cat("The next cohort will be assigned to dose level (", object$nextDose[1],",",object$nextDose[2],")", "\n")
    } else {
      if (is.na(object$overTox)) {
        cat("All tested doses are not overly toxic \n\n")
      } else {
        cat("Dose level", object$overTox, "and all levels above exhibit excessive toxicity", "\n")
      }
      if (object$decision == "stop"){
        cat("The lowest dose level is overly toxic. We terminate the entire trial for safety.")
      }else{
        cat("The current dose level is", object$curDose, "\n")
        cat("The decision regarding the direction of movement is", object$decision, "\n")
        cat("The next cohort will be assigned to dose level", object$nextDose, "\n")
      }
    }
  }
}



