#' Determination of the dose level for next cohort in the CFO-type designs with late-onset toxicities
#' 
#' Propose the next dose level in the CFO-type designs with late-onset toxicities, specifically, including 
#' Time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO) design, benchmark CFO design, TITE-aCFO design, f-aCFO 
#' design and benchmark aCFO design.
#' 
#' @usage lateonset.next(target, p.true, currdose, design, tau, enter.times, dlt.times, 
#'        current.t, doses, prior.para = list(alp.prior = target, bet.prior = 1 - target),
#'        cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05)
#'
#' @param target the target DLT rate.
#' @param p.true The true DLT rates under the different dose levels.
#' @param currdose the current dose level.
#' @param design option for selecting different designs, which can be set as \code{'TITE-CFO'}, \code{'TITE-aCFO'}, 
#'               \code{'fCFO'}, \code{'f-aCFO'}, \code{'bCFO'}, and \code{'b-aCFO'}. Specifically, \code{'bCFO'} refers 
#'               to the benchmark CFO, and \code{'b-aCFO'} denotes the benchmark aCFO.
#' @param tau maximal assessment window size
#' @param enter.times the time that each participant enters the trial.
#' @param dlt.times the time to DLT for each subject in the trial. If no DLT occurs for a subject, 
#'                  \code{dlt.times} is set to 0.
#' @param current.t the current time.
#' @param doses the dose level for each subject in the trial.
#' @param prior.para  the prior parameters for a beta distribution, usually set as list(alp.prior=target, bet.prior=1-target) by default. \code{alp.prior} 
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
#'
#' @details Late-onset outcomes commonly occur in phase I trials involving targeted agents or immunotherapies. As a 
#'          result, the TITE framework and fractional framework serve as two imputation methods to handle pending data 
#'          related to late-onset outcomes. This approach extends the original designs to integrate time information 
#'          for delayed outcomes, leading to the development of TITE-CFO, fCFO, TITE-aCFO, and f-aCFO designs. \cr
#'          In the TITE framework context, an assumption about the distribution of time to DLT must be pre-specified, 
#'          whereas the fractional framework does not require justification for a specific distribution of the time to 
#'          DLT. Consequently, fCFO and f-aCFO adapt to a more diverse range of scenarios.\cr
#'          The function \code{lateonset.next()} also provides the option to execute 
#'          the benchmark CFO and aCFO design. These two methods await complete observation of toxicity outcomes for 
#'          the previous cohorts before determining the next dose assignment. This enhances precision but comes at the 
#'          expense of a prolonged trial duration.
#' 
#' @return The \code{lateonset.next()} function returns 
#' \itemize{
#'   \item{target: }{the target DLT rate.}
#'   \item{decision: }{the decision in the CFO design, where \code{left}, \code{stay}, and \code{right} represent the 
#'   movement directions, and \code{stop} indicates stopping the experiment}
#'   \item{currdose: }{the current dose level.}
#'   \item{nextdose: }{the recommended dose level for the next cohort.}
#'   \item{overtox: }{the situation regarding which position experiences overly toxicity, where 'NA' signifies that the 
#'   occurrence of overly toxicity did not happen.}
#'   \item{tover.doses: }{a vector indicating whether the dose level (from the first to last dose level) is over-toxic 
#'   or not. }
#' }
#' 
#' @author Jialu Fang 
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Jin, H., & Yin, G. (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
#'             phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}. \cr
#'             Yin, G., Zheng, S., & Xu, J. (2013). Fractional dose-finding methods with late-onset toxicity in 
#'             phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870.
#' @import survival
#' @importFrom utils tail
#' @export
#'
#' @examples
#' ## Given the parameters for the function, the unit for time-related parameters is in months. 
#' ## The parameter generation follows the following guidelines: 
#' ## 1) the assessment window ([0, $tau]) covers three months
#' ## 2) the accrual rate is set at two patients per month. 
#' ## 3) the arrival times of patients are distributed uniformly.
#' ## 4) the time to DLT events is simulated using a Weibull distribution, with 50% of these events 
#' ##    occurring in the first half of the assessment window.
#' target <- 0.2; p.true <- c(0.02, 0.05, 0.20, 0.28, 0.34, 0.40, 0.44)
#' prior.para = list(alp.prior = target, bet.prior = 1 - target)
#' enter.times<- c(0, 0.266, 0.638, 1.54, 2.48, 3.14, 3.32, 4.01, 4.39, 5.38, 5.76,
#'                6.54, 6.66, 6.93, 7.32, 7.65, 8.14, 8.74)
#' dlt.times<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0.995, 0, 0, 0, 0, 0, 0, 0, 2.58)
#' current.t<- 9.41
#' doses<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4)
#' ## determine the dose level for the next cohort using the TITE-CFO design
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'TITE-CFO', tau = 3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' ## determine the dose level for the next cohort using the TITE-aCFO design
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'TITE-aCFO', tau = 3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' ## determine the dose level for the next cohort using the f-CFO design
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'fCFO', tau = 3, enter.times,  
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' ## determine the dose level for the next cohort using the f-aCFO design
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'f-aCFO', tau = 3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' ## determine the dose level for the next cohort using the benchmark CFO design
#' decision <- lateonset.next(target, p.true, currdose = 4, design = 'bCFO', tau = 3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' ## determine the dose level for the next cohort using the benchmark aCFO design
#' decision <- lateonset.next(target, p.true, currdose=4, design='b-aCFO', tau=3, enter.times,   
#'                dlt.times, current.t, doses, prior.para)
#' summary(decision)
#' 
lateonset.next <- function(target, p.true, currdose, design, tau, enter.times, dlt.times, 
                           current.t, doses, prior.para=list(alp.prior=target, bet.prior=1-target),
                           cutoff.eli=0.95, extrasafe=FALSE, offset=0.05){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # Below functions are to impute missing y 
  #------------------------------------------------------------------------------------------
  fracImpute <- function(enter.times, dlt.times, current.time, tau){ 
    
    #args:
    # enter.times: The enter times of the patients, a vector 
    # dlt.times: The DLT times of the patients, if no DLT, 0, a vector
    # current.time: Current time point: a value
    # tau: Observing window size
    
    #return:
    # ym: Imputed y for the patient with no DLT and follow-up time < tau
    
    assesstime = enter.times+tau;	
    dlt.times[dlt.times==0]= tau+1;
    yo = (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=(current.time-enter.times))*(current.time<assesstime);		
    No.impute <- FALSE
    if (sum(yo)==0)	{
      No.impute <- TRUE
      ym <- yo
      #stop("fraction design takes effect when at least one DLT has been observed")
    }
    if (sum(yo)!=0){			
      otime = yo*dlt.times+(1-yo)*((current.time-enter.times)*(current.time<assesstime)+tau*(assesstime<=current.time))			
      kmfit = survival::survfit(survival::Surv(otime,yo)~1)	
      ym = yo
      
      for (i in 1:length(yo)){
        if (current.time<assesstime[i] & yo[i]==0){
          ym[i]=(kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]- kmfit$surv[tail(which(kmfit$time<=tau),n=1)])/
            kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]
        }
      }
      
    }
    obsIdxs <- current.time >= assesstime
    obsIdxs[yo==1] <- TRUE
    
    
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs, No.impute=No.impute)
    res
  }
  
  TITEImpute.one <- function(followup.times, tau, y, n, prior.paras){
    #args:
    #   followup.times: The follow-up times of the pending patients at the dose level
    #   tau: Observing window size
    #   y: Num of Observed DLT at the dose level 
    #   n: Num of patients with observed results at the dose level
    #   prior.paras: a vector of 2, prior when estimating ptilde
    
    #return: 
    #  ym: imputed y 
    
    p.tilde <- (y+prior.paras[1])/(n+sum(prior.paras))
    #ym <- p.tilde * (1-followup.times/tau)
    ym <- p.tilde * (1-followup.times/tau) /((1-p.tilde)+p.tilde * (1-followup.times/tau))
    #    ym <- p.tilde * (1-followup.times/tau) /(1-p.tilde)
    #    ym[ym >1] <- 1 # trunc the value
    ym
  }
  
  TITEImpute <- function(enter.times, dlt.times, current.time, tau, dose.levels, ndose, prior.paras){
    #args:
    # enter.times: The enter times of the patients, a vector 
    # dlt.times: The DLT times of the patients, if no DLT before tau, 0, a vector
    # current.time: Current time point: a value
    # tau: Observing window size
    # dose.levels: dose level for each subject
    # ndose: num of dose levels
    # prior.paras: a vector of 2, prior when estimating ptilde
    
    #return:
    # ym: Imputed y for the patient with no DLT and follow-up time < tau
    
    assesstime = enter.times + tau;	
    followup.times <- current.time - enter.times
    dlt.times[dlt.times==0]= tau+1;
    yo <- (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=followup.times)*(current.time<assesstime);		
    obsIdxs <- current.time >= assesstime
    obsIdxs[yo==1] <- TRUE
    ym <- yo
    for (i in 1:ndose){
      doseIdxs <- dose.levels == i
      if (sum(1-obsIdxs[doseIdxs]!=0)){
        y <- sum(yo[doseIdxs])
        n <- sum(doseIdxs[obsIdxs==1])
        kpidxs <- doseIdxs & (obsIdxs!=1)
        ym.part <- TITEImpute.one(followup.times[kpidxs], tau, y, n, prior.paras)
        ym[kpidxs] <- ym.part
      }
    }
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs)
    res
  }
  
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (design == 'TITE-CFO'){accumulation = FALSE; impute.method = "TITE"
  }else if (design == 'fCFO'){accumulation = FALSE; impute.method = "frac"
  }else if (design == 'bCFO'){accumulation = FALSE; impute.method = "No"
  }else if (design == 'TITE-aCFO'){accumulation = TRUE; impute.method = "TITE"
  }else if (design == 'f-aCFO'){accumulation = TRUE; impute.method = "frac"
  }else if (design == 'b-aCFO'){accumulation = TRUE; impute.method = "No"}
  
  ndose <- length(p.true)
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  ays = NULL
  ans = NULL
  ## Obtain effective results
  if (impute.method == "frac"){
    impute.res <- fracImpute(enter.times, dlt.times, current.t, tau)
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    if (impute.res$No.impute){
      for (i in 1:ndose){
        ays <- c(ays, sum(y.raw[doses==i]))
        ans <- c(ans, sum(doses==i))
      } 
    }
    else{
      for (i in 1:ndose){
        ays <- c(ays, sum(y.impute[doses==i]))
        ans <- c(ans, sum(doses==i))
      }
    }
  }else if(impute.method == "TITE"){
    impute.res <-  TITEImpute(enter.times, dlt.times, current.t, tau, doses, ndose, c(target/2, 1-target/2))
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    for (i in 1:ndose){
      ays <- c(ays, sum(y.impute[doses==i]))
      ans <- c(ans, sum(doses==i))
    }
  }else if(impute.method == "No"){
    assesstime = enter.times+tau;	
    dlt.times[dlt.times==0]= tau+1;
    y.impute <- (dlt.times<=tau)*(assesstime<=current.t)
    for (i in 1:ndose){
      ays <- c(ays, sum(y.impute[doses==i]))
      ans <- c(ans, sum(doses==i))
    }
  }
  
  tover.doses <- rep(0, ndose)
  
  
  for (i in 1:ndose){
    cy <- ays[i]
    cn <- ans[i]
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[i:ndose] <- 1
      break()
    }
  }

  if (extrasafe) {
    cy <- ays[1]
    cn <- ans[1]
    prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, cutoff.eli-offset, prior.para)){
      tover.doses[1:ndose] <- 1
    }
  }
  
  position <- which(tover.doses == 1)[1]
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if (accumulation == FALSE){
    if (currdose==1){
      cys <- c(NA, ays[1:(currdose+1)])
      cns <- c(NA, ans[1:(currdose+1)])
    }else if (currdose==ndose){
      cys <- c(ays[(currdose-1):ndose], NA)
      cns <- c(ans[(currdose-1):ndose], NA)
    }else {
      cys <- ays[(currdose-1):(currdose+1)]
      cns <- ans[(currdose-1):(currdose+1)]
    }
    res <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, extrasafe, offset)
  }else{
    res <- aCFO.next (target, ays, ans, currdose, prior.para, cutoff.eli, extrasafe, offset)
  }
  nextdose <- res$nextdose
  decision <- res$decision
  overtox <- res$overtox

  out <- list(target=target, ays=ays, ans=ans, decision=decision, currdose = currdose, 
              nextdose=nextdose, overtox=overtox, tover.doses=tover.doses)
  class(out) <- "cfo"
  return(out)
}


