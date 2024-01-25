#' Conduct one simulation using the Calibration-Free Odds (CFO) or accumulative CFO (aCFO) design and find the maximum tolerated dose (MTD).
#' 
#' Use this function to conduct one simulation using the Calibration-Free Odds (CFO) or accumulative CFO (aCFO) design and find the maximum tolerated dose (MTD).
#'
#' @usage CFO.simu(target, p.true, ncohort, init.level=1, cohortsize=3, design,
#'        prior.para=list(alp.prior=target, bet.prior=1-target), seed=NULL,
#'        cutoff.eli=0.95, extrasafe=FALSE, offset=0.05)
#'
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param ncohort the total number of cohorts.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param cohortsize the number of patients or size of each cohort. The default value \code{cohortsize} is 3.
#' @param design option for selecting different designs, which can be set as \code{'CFO'} and \code{'aCFO'}.
#' @param prior.para the prior parameters for a beta distribution, usually set as list(alp.prior=target, bet.prior=1-target) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe set \code{extrasafe=TRUE} to impose a more strict early stopping rule for
#'                   extra safety
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'                stopping rule is when \code{extrasafe=TRUE}. A larger value leads to
#'                a more strict stopping rule. The default value \code{offset=0.05}
#'                generally works well.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results, the default is set to NULL.
#'                            
#' @details The \code{CFO.simu()} function is designed to determine the maximum tolerated dose (MTD) for a single CFO or aCFO 
#'          trial. If \code{design = 'CFO'}, this trial corresponds to the CFO design. If \code{design = 'aCFO'}, it
#'          corresponds to the aCFO design. \cr
#'          Given the toxicity outcomes from previous cohorts, each cohort is sequentially assigned to the most suitable dose 
#'          level based on the CFO or aCFO decision rule. An early stopping and dose elimination are incorporated into the CFO or aCFO design 
#'          to ensure patient safety and benefit. If there is substantial evidence indicating that the current dose level 
#'          exhibits excessive toxicity (\eqn{\Pr(p_C > \target|x_C, m_C \geq 3) > 0.95}), we exclude the current dose level as 
#'          well as higher dose levels from the trial. If the lowest dose level is overly toxic, the trial will be terminated 
#'          according to the early stopping rule. Upon the predefined maximum sample size is reached or the lowest dose 
#'          level is overly toxicity, the experiment is concluded, and the MTD is determined using isotonic regression.
#'
#' @return The \code{CFO.simu} function returns a list object comprising the following components:
#'         \itemize{
#'         \item{target: }{the target DLT rate.}
#'         \item{MTD: }{the selected MTD. \code{MTD=99} indicates that this trial is terminated due to early stopping.}
#'         \item{correct: }{a binary indicator of whether the recommended dose level matches the target DLT rate (1 for yes).}
#'         \item{npatients: }{the total number of patients allocated for all dose levels}
#'         \item{ntox: }{the total number of DLTs observed for all dose levels.}
#'         \item{npercent: }{the percentage of subjects assigned to the target DLT rate.}
#'         \item{over.doses: }{a vector indicating whether each dose is overdosed or not (1 for yes).}
#'         \item{cohortdose: }{a vector including the dose level assigned for each cohort.}
#'         \item{patientDLT: }{a vector including the DLT outcome observed for each patient.}
#'         }
#' 
#' @author Jialu Fang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#'
#' @examples
#' target <- 0.2; ncohort <- 12; cohortsize <- 3
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' ## find the MTD for a single CFO trial
#' CFOtrial <- CFO.simu(target, p.true, ncohort, init.level=1, cohortsize=3, design='CFO',
#'             prior.para=list(alp.prior=target, bet.prior=1-target))
#' summary(CFOtrial)
#' plot(CFOtrial)
#' ## find the MTD for a single aCFO trial
#' aCFOtrial <- CFO.simu(target, p.true, ncohort, init.level=1, cohortsize=3, design='aCFO',
#'              prior.para=list(alp.prior=target, bet.prior=1-target))
#' summary(aCFOtrial)
#' plot(aCFOtrial)
#' @export
CFO.simu <- function(target, p.true, ncohort, init.level=1, cohortsize=3, design,
                     prior.para=list(alp.prior=target, bet.prior=1-target), seed=NULL,
                     cutoff.eli=0.95, extrasafe=FALSE, offset=0.05){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
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
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  set.seed(seed)
  earlystop <- 0
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  currdose <- init.level
  
  ays <- rep(0, ndose) # number of responses for different doses.
  ans <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
  DLTlist <- c()
  
  for (i in 1:ncohort){
    pc <- p.true[currdose]
    doselist[i] <- currdose
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    DLTlist <- c(DLTlist, cres)
    
    # update results
    ays[currdose] <- ays[currdose] + sum(cres)
    ans[currdose] <- ans[currdose] + cohortsize
    
    cy <- ays[currdose]
    cn <- ans[currdose]
    
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[currdose:ndose] <- 1
    }
    
    if (currdose == 1){
      if (extrasafe) {
        cy <- ays[1]
        cn <- ans[1]
        prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
        if (overdose.fn(target, cutoff.eli, prior.para)){
          tover.doses[1:ndose] <- 1
        }
      }
    }
    
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    }
    
    prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (design == 'CFO'){
      # the results for current 3 dose levels
      if (currdose!=1){
        cys <- ays[(currdose-1):(currdose+1)]
        cns <- ans[(currdose-1):(currdose+1)]
      }else{
        cys <- c(NA, ays[1:(currdose+1)])
        cns <- c(NA, ans[1:(currdose+1)])
      }
      currdose <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, extrasafe, offset)$nextdose
    }else if (design == 'aCFO'){
      currdose <- aCFO.next(target, ays, ans, currdose, prior.para, cutoff.eli, extrasafe, offset)$nextdose
    }else{
      stop("The input design is invalid; it can only be set as 'CFO' or 'aCFO'.")
    }
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
  
  out<-list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, npercent=npercent, 
            over.doses=tover.doses, cohortdose=doselist, patientDLT=DLTlist)
  class(out) <- "cfo"
  return(out)
}
