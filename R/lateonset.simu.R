#' Find the maximum tolerated dose (MTD) for a single trial using the CFO-type designs with late-onset toxicities.
#'
#' Use this function to find the maximum tolerated dose (MTD) for the CFO-type designs with late-onset toxicities, 
#' specifically, including Time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO), benchmark CFO, TITE-aCFO design, 
#' the f-aCFO design and benchmark aCFO design.
#'
#' @usage lateonset.simu(target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'        accrual.dist, design, init.level = 1, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        seed = NULL, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05)
#'
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param tau maximal assessment window size.
#' @param cohortsize the number of patients or size of each cohort. The default value \code{cohortsize} is 3.
#' @param ncohort the total number of cohorts.
#' @param tte.para the parameter related with the distribution of the time to DLT events. The time to DLT is sample from a Weibull 
#'                 distribution, with \code{tte.para} representing the proportion of DLTs occurring within the initial half of the 
#'                 assessment window \code{tau}.
#' @param accrual the accrual rate, i.e., the number of patients accrued in per unit time.
#' @param accrual.dist the distribution of the arrival times of patients. When \code{accrual.dist = 'fix'}, it corresponds to all 
#'                     patients in each cohort arriving simultaneously at a given accrual rate. When \code{accrual.dist = 'unif'}, 
#'                     it corresponds to a uniform distribution, and when \code{accrual.dist = 'exp'}, it corresponds to an 
#'                     exponential distribution.
#' @param design option for selecting different designs, which can be set as \code{'TITE-CFO'}, \code{'TITE-aCFO'}, 
#'               \code{'fCFO'}, \code{'f-aCFO'}, \code{'bCFO'}, and \code{'b-aCFO'}. Specifically, \code{'bCFO'} refers 
#'               to the benchmark CFO, and \code{'b-aCFO'} denotes the benchmark aCFO.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param prior.para  the prior parameters for a beta distribution, usually set as list(alp.prior = target, bet.prior = 1 - target) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli = 0.95}) for general use.
#' @param extrasafe set \code{extrasafe = TRUE} to impose a more strict early stopping rule for
#'                   extra safety.
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'                stopping rule is when \code{extrasafe=TRUE}. A larger value leads to
#'                a more strict stopping rule. The default value \code{offset = 0.05}
#'                generally works well.
#' @param seed an integer to set as the seed of the random number generator for reproducible results, the default is set to NULL.
#' @note  The results returned by \code{lateonset.simu()} represent a single trial outcome. In practical 
#'        applications, multiple simulations are required to ascertain the MTD, and this can be achieved using the 
#'        \code{CFO.oc()} function.
#'
#' @return The \code{lateonset.simu()} function returns a list object comprising the following components: 
#' \itemize{
#' \item{target: }{the target DLT rate.}
#' \item{MTD: }{the selected MTD. \code{MTD=99} indicates that this trial is terminated due to early stopping.}
#' \item{correct: }{a binary indicator of whether the recommended dose level matches the target DLT rate (1 for yes).}
#' \item{npatients: }{the total number of patients allocated for all dose levels}
#' \item{ntox: }{the total number of DLTs observed for all dose levels.}
#' \item{npercent: }{the percentage of subjects assigned to the target DLT rate.}
#' \item{over.doses: }{a vector indicating whether each dose is overdosed or not (1 for yes).}
#' \item{cohortdose: }{a vector including the dose level assigned for each cohort.}
#' \item{patientDLT: }{a vector including the DLT outcome observed for each patient.}
#' \item{totaltime: }{the duration of the trial.}
#' \item{entertimes: }{the time that each participant enters the trial.}
#' \item{DLT.times: }{the time to DLT for each subject in the trial. If no DLT occurs for a certain subject, 
#'                  \code{DLT.times} is 0.}
#' }
#' 
#'         
#' @author Jialu Fang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Jin, H., & Yin, G. (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
#'             phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}. \cr
#'             Yin, G., Zheng, S., & Xu, J. (2013). Fractional dose-finding methods with late-onset toxicity in 
#'             phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870.
#' @export
#'
#' @examples
#' target <- 0.2; ncohort <- 12; cohortsize <- 3
#' p.true <- c(0.01, 0.02, 0.05, 0.20, 0.30, 0.50, 0.70)
#' tau <- 3; accrual <- 2; tte.para <- 0.5; accrual.dist <- 'unif'
#' ## find the MTD for a single TITE-CFO trial
#' TITECFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'TITE-CFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(TITECFOtrial)
#' plot(TITECFOtrial)
#' ## find the MTD for a single TITE-aCFO trial
#' TITEaCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'TITE-aCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(TITEaCFOtrial)
#' plot(TITEaCFOtrial)
#' ## find the MTD for a single fCFO trial
#' fCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'fCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(fCFOtrial)
#' plot(fCFOtrial)
#' ## find the MTD for a single f-aCFO trial
#' faCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'f-aCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(faCFOtrial)
#' plot(faCFOtrial)
#' ## find the MTD for a single benchmark CFO trial
#' bCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'bCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(bCFOtrial)
#' plot(bCFOtrial)
#' ## find the MTD for a single benchmark aCFO trial
#' baCFOtrial <- lateonset.simu (target, p.true, tau, cohortsize, ncohort, tte.para, accrual, 
#'                 accrual.dist, design = 'b-aCFO', init.level = 1, 
#'                 prior.para = list(alp.prior = target, bet.prior = 1 - target), seed = 1)
#' summary(baCFOtrial)
#' plot(baCFOtrial)
lateonset.simu <- function(target, p.true, tau, cohortsize, ncohort, tte.para, accrual, accrual.dist, 
                    design, init.level=1, prior.para=list(alp.prior=target, bet.prior=1-target), seed=NULL,
                    cutoff.eli=0.95, extrasafe=FALSE, offset=0.05){
  
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  
  # The function is to obtain the DLT results (with TITE) for each subject
  gen.tite<-function(n, pi, tau=1, alpha=0.5){
    #args:
    #   n: Num of subjects to generate
    #   pi: Target DLT rate, pi=Pr(T<=tau)
    #   tau: Maximal window size
    #   alpha: Parameter for generate time
    #Return:
    #   if no DLT, tox.t=0
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1-pi)/log(1-pihalft))/log(2);
      lambda = -log(1-pi)/(tau^alpha);
      t = (-log(runif(n))/lambda)^(1/alpha);
      return(t);
    }
    ############ end of subroutines ############
    tox = rep(0, n);
    t.tox = rep(0, n);
    #### Weibull
    pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = weib(n, pi, pihalft);
    tox[t.tox<=tau]=1;
    ntox.st = sum(tox);
    t.tox[tox==0]=0;
    return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ############################################################################### 
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  set.seed(seed)
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  
  earlystop <- 0
  enter.times <- NULL # enter time of each subject
  dlt.times <- NULL # dlt time of each subject
  dlts <- NULL # dlt event for each subject
  doses <- NULL # dose level for each subject
  current.t<- 0
  currdose <- init.level  #current dose level
  
  tover.doses <- rep(0, ndose)
  
  for (i in 1:ncohort){
    curP <- p.true[currdose]
    doselist[i] <- currdose
    
    if (accrual.dist=='fix'){
      delta.times <- rep(0, cohortsize)
    }else if (accrual.dist == 'unif'){
      delta.times <- cumsum(c(0, runif(cohortsize-1, 0, 2/accrual)))
    }else if (accrual.dist == 'exp'){
      delta.times <- cumsum(c(0, rexp(cohortsize-1, rate=accrual)))
    }

    enter.times <- c(enter.times, current.t+delta.times)
    
    # obtain the results of the patients
    obscohort <- gen.tite(cohortsize, curP, alpha=tte.para, tau=tau);
    dlt.times <- c(dlt.times, obscohort$t.tox);
    dlts <- c(dlts, obscohort$tox);
    doses <- c(doses, rep(currdose, cohortsize));
    
    # Move to next cohort 
    if (i != ncohort){
      if (accrual.dist=='fix'){
        delta.time <- cohortsize/accrual
      }else if (accrual.dist == 'unif'){
        delta.time <- runif(1, 0, 2/accrual)
      }else if (accrual.dist == 'exp'){
        delta.time <- rexp(1, rate=accrual)
      }
    }else{
      delta.time <- tau
    }
    current.t<- enter.times[length(enter.times)] + delta.time;
    
    if (design == 'bCFO' || design == 'b-aCFO'){
      current.t <- enter.times[length(enter.times)] + tau
      res <- lateonset.next(target, p.true, currdose, design, tau, enter.times, dlt.times, current.t, doses, 
                            prior.para, cutoff.eli, extrasafe, offset)
      tover.doses <- res$tover.doses
      overTox <- res$overTox
      current.t <- current.t + delta.time
    }else{
      res <- lateonset.next(target, p.true, currdose, design, tau, enter.times, dlt.times, current.t, doses, 
                            prior.para, cutoff.eli, extrasafe, offset)
      tover.doses <- res$tover.doses
      overTox <- res$overTox
    }
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    } else{
      currdose <- res$nextdose
    }
  }
  
  ans <- NULL
  ays <- NULL
  assess.t <- enter.times + tau
  y.raw <- (dlt.times!=0)*1
  for (j in 1:ndose){
    ans <- c(ans, sum(doses==j))
    ays <- c(ays, sum(y.raw[doses==j]))
  }
  if (earlystop==0){
    MTD <- select.mtd(target, ans, ays, prior.para, cutoff.eli, extrasafe, offset, verbose=FALSE)$MTD
  }else{
    MTD <- 99
  }
  
  correct <- 0
  if (MTD == target){
    correct <- 1
  }
  
  npercent <- ans[which(p.true == target)]/(ncohort*cohortsize)
  
  out <- list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, npercent=npercent, 
              over.doses=tover.doses, cohortdose=doselist, patientDLT = dlts, 
              totaltime=assess.t[length(assess.t)], entertimes=enter.times, DLTtimes=dlt.times)
  class(out) <- "cfo"
  return(out)
}
