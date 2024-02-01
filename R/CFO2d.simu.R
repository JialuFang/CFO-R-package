#' Obtain the operating characteristics of the 2dCFO designs for multiple simulations.
#'
#' @usage CFO2d.simu(target, p.true, init.level=c(1,1), ncohort, cohortsize, 
#'                  prior.para=list(alp.prior=target, bet.prior=1-target),
#'                  cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, seed=NULL)
#'
#' @param target the target DLT rate.
#' @param p.true a matrix representing the true DIL rates under the different dose levels.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is c(1,1).
#' @param ncohort the total number of cohorts.
#' @param cohortsize the number of patients or size of each cohort. 
#' @param prior.para the prior parameters for a beta distribution, usually set as list(alp.prior=target, bet.prior=1-target) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior})
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli = 0.95}) for general use.
#' @param extrasafe set \code{extrasafe = TRUE} to impose a more strict early stopping rule for
#'                   extra safety.
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'                stopping rule is when \code{extrasafe=TRUE}. A larger value leads to
#'                a more strict stopping rule. The default value \code{offset = 0.05}
#'                generally works well.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results, the default is set to NULL.
#'
#' @details 
#' The `CFO2d.simu` function simulates the operating characteristics of the CFO designs 
#' in a two-dimensional dose-finding trial. The function uses parameters such as target DLT rate, true DLT rates 
#' under different dose levels, and cohort details, and provides detailed output for performance evaluation. 
#' @author Wenliang Wang
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item{target: }{the target DLT rate.}
#'   \item{MTD: }{a vector of length 2 representing the recommended dose level. \code{MTD=99} indicates that this trial is terminated due to early stopping.}
#'   \item{correct: }{a binary indicator of whether the recommended dose level matches the target DLT rate (1 for yes).}
#'   \item{npatients: }{a matrix of the number of patients allocated for different doses.}
#'   \item{ntox: }{a matrix of the number of DLT observed for different doses.}
#'   \item{npercent: }{the percentage of patients assigned to the target DLT rate.}
#'   \item{over.doses: }{a matrix indicating whether each dose is overdosed or not (1 for yes).}
#'   \item{cohortdose: }{the dose combination assigned for each cohort.}
#'   \item{ptoxic: }{The percentage of subjects assigned to dose levels with a DLT rate greater than the target.}
#'   \item{patientDLT: }{the DLT observed at each cohort.}
#'   \item{sumDLT: }{the total number of DLT observed.}
#'   \item{earlystop: }{a binary indicator of whether the trial is early stopped (1 for yes).}
#' }
#' @export
#'
#' @examples
#' ## Simulate a two-dimensional dose-finding trial with 20 cohorts of size 3.
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#' target <- 0.3; ncohort <- 20; cohortsize <- 3
#' CFO2dtrial <- CFO2d.simu(target, p.true, init.level = c(1,1), ncohort, cohortsize, seed = 1)
#' summary(CFO2dtrial)
#' plot(CFO2dtrial)



CFO2d.simu <- function(target, p.true, init.level=c(1,1), ncohort, cohortsize,  
                       prior.para=list(alp.prior=target, bet.prior=1-target), 
                       cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, seed=NULL){
  # target: Target DIL rate
  # p.true: True DIL rates under the different dose levels
  # ncohort: The number of cohorts
  # cohortsize: The sample size in each cohort
  # alp.prior, bet.prior: prior parameters
  if (!is.null(seed)){
    set.seed(seed)
  }
  earlystop <- 0
  ndose.A <- length(p.true[,1])
  ndose.B <- length(p.true[1,])
  cidx.A <- init.level[1]
  cidx.B <- init.level[2]
  obs <- list()
  ays <- matrix(0, ndose.A, ndose.B) # number of responses for different doses.
  ans <- matrix(0, ndose.A, ndose.B) # number of subject for different doses.
  tover.doses <- matrix(0, ndose.A, ndose.B) # Whether each dose is overdosed or not, 1 yes
  
  # Initialize vectors to store dose combinations and number of DLTs for each cohort
  # simu.res.dose <- vector("list", ncohort) # Change to list to store dose pairs
  simu.res.dose <- matrix(nrow = ncohort, ncol = 2)
  simu.res.DLT <- matrix(nrow = ncohort, ncol = cohortsize)
  
  overdose.2d <- function(phi, threshold, obs, prior.para=list(alp.prior=phi, bet.prior=1-phi)){
    y <- obs$y
    n <- obs$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    if ((pp >= threshold) & (obs$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  post.prob.fn <- function(phi, y, n, alp.prior=phi, bet.prior=1-phi){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  for (i in 1:ncohort){
    pc <- p.true[cidx.A, cidx.B] 
    if (!is.null(seed)) {
      iter_seed <- (seed * 100) + i
      set.seed(iter_seed)
    }
    cres <- rbinom(cohortsize, 1, pc)
    ays[cidx.A, cidx.B] <- ays[cidx.A, cidx.B] + sum(cres)
    ans[cidx.A, cidx.B] <- ans[cidx.A, cidx.B] + cohortsize
    
    simu.res.dose[i, ] <- c(cidx.A, cidx.B)
    simu.res.DLT[i,] <- cres
    
    cy <- ays[cidx.A, cidx.B]
    cn <- ans[cidx.A, cidx.B]
    
    obs <- c(list(y=cy, n=cn, ays=ays, ans=ans, cidx.A=cidx.A, cidx.B=cidx.B), obs)
    
    if (overdose.2d(target, cutoff.eli, obs)){
      tover.doses[cidx.A:ndose.A, cidx.B:ndose.B] <- 1
    }
    if (cidx.A == 1 & cidx.B == 1) {
      if (extrasafe) {
        if (overdose.2d(target, cutoff.eli-offset, obs)){
          tover.doses[1:1] <- 1
        }
      }
    }
    if (tover.doses[1,1] == 1){
      earlystop <- 1
      break()
    }
    
    if (cidx.A!=1 & cidx.B!=1 & cidx.A!=ndose.A & cidx.B!=ndose.B){
      # no boundary
      cys <- ays[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cns <- ans[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cover.doses <- tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
    } else if (cidx.A==1 & cidx.B==1){
      # (1, 1)
      cys <- rbind(c(NA,NA,NA),cbind(c(NA,NA),ays[1:2,1:2]))
      cns <- rbind(c(NA,NA,NA),cbind(c(NA,NA),ans[1:2,1:2]))
      cover.doses <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tover.doses[1:2,1:2]))
    } else if (cidx.A==ndose.A & cidx.B==ndose.B){
      # (nA, nB)
      cys <- rbind(cbind(ays[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cns <- rbind(cbind(ans[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cover.doses <- rbind(cbind(tover.doses[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B==ndose.B){
      # (1, nB) 
      cys <- rbind(c(NA,NA,NA),cbind(ays[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cns <- rbind(c(NA,NA,NA),cbind(ans[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cover.doses <- rbind(c(NA,NA,NA),cbind(tover.doses[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
    } else if (cidx.A==ndose.A & cidx.B==1){
      # (nA, 1) 
      cys <- rbind(cbind(c(NA,NA), ays[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cns <- rbind(cbind(c(NA,NA), ans[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cover.doses <- rbind(cbind(c(NA,NA), tover.doses[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B!=1){
      # (1, 2:(nB-1))
      cys <- rbind(c(NA,NA,NA), ays[1:2, (cidx.B-1):(cidx.B+1)])
      cns <- rbind(c(NA,NA,NA), ans[1:2, (cidx.B-1):(cidx.B+1)])
      cover.doses <- rbind(c(NA,NA,NA), tover.doses[1:2, (cidx.B-1):(cidx.B+1)])
    } else if (cidx.A!=1 & cidx.B==1){
      # (2:(nA-1), 1)
      cys <- cbind(c(NA,NA,NA), ays[(cidx.A-1):(cidx.A+1), 1:2])
      cns <- cbind(c(NA,NA,NA), ans[(cidx.A-1):(cidx.A+1), 1:2])
      cover.doses <- cbind(c(NA,NA,NA), tover.doses[(cidx.A-1):(cidx.A+1), 1:2])
    } else if (cidx.A==ndose.A & cidx.B!=ndose.B){
      # (nA, 2:(nB-1))
      cys <- rbind(ays[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cns <- rbind(ans[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cover.doses <- rbind(tover.doses[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
    } else if (cidx.A!=ndose.A & cidx.B==ndose.B){
      # (2:(nA-1), nB)
      cys <- cbind(ays[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cns <- cbind(ans[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cover.doses <- cbind(tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
    } else {
      message('no such case')
    }
    
    idx.chg <- CFO2d.next(target, cys, cns, c(cidx.A, cidx.B), prior.para, seed=seed)$index
    cidx.A <- cidx.A + idx.chg[1]
    cidx.B <- cidx.B + idx.chg[2]
  }
  if (earlystop==0){
    MTD <- CFO2d.selectmtd(target, ans, ays)$MTD
  }else{
    MTD <- c(99,99)
  }
  
  correct <- 0
  if(MTD[1]>ndose.A | MTD[2]>ndose.B){
    correct <- 0
  } else if (length(MTD)!=2){
    correct <- 0
  }else if (p.true[MTD[1],MTD[2]]==target){
    correct <- 1
  }
  
  npercent <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]==target){
        npercent <- npercent + ans[j,k]
      }
    }
  }
  npercent <- npercent/(ncohort*cohortsize)
  
  ptoxic <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]>target){
        ptoxic <- ptoxic + ans[j,k]
      }
    }
  }
  ptoxic <- ptoxic/(ncohort*cohortsize)
  # simu.res <- list(dose = simu.res.dose, DLT = simu.res.DLT)
  out<-list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, 
            npercent=npercent, over.doses=tover.doses, cohortdose=simu.res.dose, ptoxic=ptoxic,
            patientDLT = simu.res.DLT, sumDLT=sum(simu.res.DLT), earlystop=earlystop)
  class(out) <- "cfo"
  return(out)
}
