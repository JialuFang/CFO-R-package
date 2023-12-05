#' Generate operating characteristics for multiple simulations
#' 
#' Obtain the operating characteristics of the CFO-type and aCFO-type designs for multiple simulations.
#'
#' @usage CFO2d.sim(phi, p.true, ncohort=20, cohortsize=3, init.level=c(1,1), 
#'                  add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL)
#'
#' @param phi the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param ncohort the total number of cohorts, the default value is 20.
#' @param cohortsize the sample size in each cohort, the default value is 3. 
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is c(1,1).
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param seed an integer to set as the seed of the random number generator for reproducible results, the default is set to NULL.
#'
#' @details 
#' The `CFO2d.sim` function simulates the operating characteristics of the CFO designs 
#' in a two-dimensional dose-finding trial. The function uses parameters such as target DLT rate, true DLT rates 
#' under different dose levels, and cohort details, and provides detailed output for performance evaluation. 
#' It relies on the BOIN package.
#' @author Wenliang Wang
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item{MTD: }{A vector of length 2 representing the recommended dose level.}
#'   \item{dose.ns: }{A matrix of the number of patients allocated for different doses.}
#'   \item{DLT.ns: }{A matrix of the number of DLT observed for different doses.}
#'   \item{p.true: }{The matrix of the true DLT rates under the different dose levels.}
#'   \item{target: }{The target DLT rate.}
#'   \item{over.doses: }{A matrix indicating whether each dose is overdosed or not (1 for yes).}
#'   \item{correct: }{A binary indicator of whether the recommended dose level matches the target DLT rate.}
#'   \item{npercent: }{The percentage of subjects assigned to the target DLT rate.}
#'   \item{ptoxic: }{The percentage of subjects assigned to dose levels with a DLT rate greater than the target.}
#'   \item{ntox: }{The total number of DLTs observed.}
#'   \item{dose: }{The dose combination assigned for each cohort.}
#'   \item{DLT: }{The DLT observed at each cohort.}
#' }
#' @import BOIN
#' @export
#'
#' @examples
#' ## Simulate a two-dimensional dose-finding trial with 20 cohorts of size 3.
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#'
#' CFO2d.sim(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3)



CFO2d.sim <- function(phi, p.true, ncohort=20, cohortsize=3, init.level=c(1,1), add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL){
  # phi: Target DIL rate
  # p.true: True DIL rates under the different dose levels
  # ncohort: The number of cohorts
  # cohortsize: The sample size in each cohort
  # alp.prior, bet.prior: prior parameters
  set.seed(seed)
  earlystop <- 0
  ndose.A <- length(p.true[,1])
  ndose.B <- length(p.true[1,])
  cidx.A <- init.level[1]
  cidx.B <- init.level[2]
  obs <- list()
  tys <- matrix(0, ndose.A, ndose.B) # number of responses for different doses.
  tns <- matrix(0, ndose.A, ndose.B) # number of subject for different doses.
  tover.doses <- matrix(0, ndose.A, ndose.B) # Whether each dose is overdosed or not, 1 yes
  
  # Initialize vectors to store dose combinations and number of DLTs for each cohort
  # sim.res.dose <- vector("list", ncohort) # Change to list to store dose pairs
  sim.res.dose <- matrix(nrow = ncohort, ncol = 2)
  sim.res.DLT <- vector("numeric", ncohort)
  
  overdose.fn <- function(phi, obs, add.args=list(alp.prior=phi, bet.prior=1-phi)){
    y <- obs$y
    n <- obs$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    if ((pp >= 0.95) & (obs$n>=3)){
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
    cres <- rbinom(cohortsize, 1, pc)
    tys[cidx.A, cidx.B] <- tys[cidx.A, cidx.B] + sum(cres)
    tns[cidx.A, cidx.B] <- tns[cidx.A, cidx.B] + cohortsize
    
    # sim.res.dose[[i]] <- c(cidx.A, cidx.B) # Store as a pair
    sim.res.dose[i, ] <- c(cidx.A, cidx.B)
    sim.res.DLT[i] <- sum(cres)
    
    cy <- tys[cidx.A, cidx.B]
    cn <- tns[cidx.A, cidx.B]
    
    obs <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx.A=cidx.A, cidx.B=cidx.B), obs)
    
    if (overdose.fn(phi, obs)){
      tover.doses[cidx.A:ndose.A, cidx.B:ndose.B] <- 1
    }
    
    if (tover.doses[1,1] == 1){
      earlystop <- 1
      break()
    }
    
    if (cidx.A!=1 & cidx.B!=1 & cidx.A!=ndose.A & cidx.B!=ndose.B){
      # no boundary
      cys <- tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cns <- tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
      cover.doses <- tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
    } else if (cidx.A==1 & cidx.B==1){
      # (1, 1)
      cys <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tys[1:2,1:2]))
      cns <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tns[1:2,1:2]))
      cover.doses <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tover.doses[1:2,1:2]))
    } else if (cidx.A==ndose.A & cidx.B==ndose.B){
      # (nA, nB)
      cys <- rbind(cbind(tys[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cns <- rbind(cbind(tns[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
      cover.doses <- rbind(cbind(tover.doses[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B==ndose.B){
      # (1, nB) 
      cys <- rbind(c(NA,NA,NA),cbind(tys[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cns <- rbind(c(NA,NA,NA),cbind(tns[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
      cover.doses <- rbind(c(NA,NA,NA),cbind(tover.doses[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
    } else if (cidx.A==ndose.A & cidx.B==1){
      # (nA, 1) 
      cys <- rbind(cbind(c(NA,NA), tys[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cns <- rbind(cbind(c(NA,NA), tns[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
      cover.doses <- rbind(cbind(c(NA,NA), tover.doses[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
    } else if (cidx.A==1 & cidx.B!=1){
      # (1, 2:(nB-1))
      cys <- rbind(c(NA,NA,NA), tys[1:2, (cidx.B-1):(cidx.B+1)])
      cns <- rbind(c(NA,NA,NA), tns[1:2, (cidx.B-1):(cidx.B+1)])
      cover.doses <- rbind(c(NA,NA,NA), tover.doses[1:2, (cidx.B-1):(cidx.B+1)])
    } else if (cidx.A!=1 & cidx.B==1){
      # (2:(nA-1), 1)
      cys <- cbind(c(NA,NA,NA), tys[(cidx.A-1):(cidx.A+1), 1:2])
      cns <- cbind(c(NA,NA,NA), tns[(cidx.A-1):(cidx.A+1), 1:2])
      cover.doses <- cbind(c(NA,NA,NA), tover.doses[(cidx.A-1):(cidx.A+1), 1:2])
    } else if (cidx.A==ndose.A & cidx.B!=ndose.B){
      # (nA, 2:(nB-1))
      cys <- rbind(tys[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cns <- rbind(tns[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
      cover.doses <- rbind(tover.doses[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
    } else if (cidx.A!=ndose.A & cidx.B==ndose.B){
      # (2:(nA-1), nB)
      cys <- cbind(tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cns <- cbind(tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
      cover.doses <- cbind(tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
    } else {
      message('no such case')
    }
    
    idx.chg <- CFO2d.next(phi, cys, cns, cover.doses, c(cidx.A, cidx.B), add.args, seed)$decision
    cidx.A <- cidx.A + idx.chg[1]
    cidx.B <- cidx.B + idx.chg[2]
  }
  
  est_p <- select.mtd.comb(phi, tns, tys, boundMTD=TRUE)
  if (earlystop==0){
    MTD <- select.mtd.comb(phi, tns, tys)$MTD
  }else{
    MTD <- c(99,99)
  }
  
  correct <- 0
  if(MTD[1]>ndose.A | MTD[2]>ndose.B){
    correct <- 0
  } else if (length(MTD)!=2){
    correct <- 0
  }else if (p.true[MTD[1],MTD[2]]==phi){
    correct <- 1
  }
  
  npercent <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]==phi){
        npercent <- npercent + tns[j,k]
      }
    }
  }
  npercent <- npercent/(ncohort*cohortsize)
  
  ptoxic <- 0
  for (j in 1:ndose.A) {
    for (k in 1:ndose.B) {
      if (p.true[j,k]>phi){
        ptoxic <- ptoxic + tns[j,k]
      }
    }
  }
  ptoxic <- ptoxic/(ncohort*cohortsize)
  # sim.res <- list(dose = sim.res.dose, DLT = sim.res.DLT)
  list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, over.doses=tover.doses, correct=correct, npercent=npercent, ptoxic=ptoxic, ntox=sum(tys), dose=sim.res.dose, DLT = sim.res.DLT)
}


# p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
# 0.10, 0.15, 0.30, 0.45, 0.55,
# 0.15, 0.30, 0.45, 0.50, 0.60),
# nrow = 3, ncol = 5, byrow = TRUE)
# 
# CFO2d.res <- CFO2d.sim(phi=0.3, p.true, ncohort = 20, cohortsize = 3)


