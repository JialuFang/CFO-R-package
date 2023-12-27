#' 
#' Determination of the dose level for next cohort
#' 
#' Use the function to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
#'
#' @usage CFO.next(phi, cys, cns, curDose, add.args=list(alp.prior=phi,bet.prior=1-phi), seed=NULL)
#'
#' @param phi the target DLT rate
#' @param cys the current number of DLTs observed at the left, current, and right dose levels.
#' @param cns the current number of patients treated at the left, current, and right dose levels.
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param curDose the current dose level.
#' @param seed an integer to set as the seed of the random number generator for reproducible results.
#'
#' @details The CFO design determines the dose level for the next cohort by assessing evidence from the current 
#'          dose level and its adjacent levels. This evaluation is based on odds ratios denoted as \eqn{O_k}, where 
#'          k = L, C, R represents left, current, and right dose levels. Additionally, we define \eqn{\overline{O}_k = 1/O_k}. 
#'          The ratio \eqn{O_C / \overline{O}_{L}} indicates the inclination for de-escalation, while \eqn{\overline{O}_C / O_R} 
#'          quantifies the tendency for escalation. Threshold values \eqn{\gamma_L} and \eqn{\gamma_R} are chosen to 
#'          minimize the probability of making incorrect decisions.The decision process is summarized in Table 1
#'          of Jin and Yin (2022).
#'          An overdose control rule is implemented to ensure patient safety. If the data suggest excessive 
#'          toxicity at the current dose level, we exclude that level and those higher levels. Two scenarios 
#'          lead to a decision on one side only: when the current dose is at the boundary (the first or last dose level) 
#'          or when higher dose levels have been eliminated.
#'          
#' @note    When the current dose level is the lowest or highest (i.e., at the boundary), the parts in \code{cys} and 
#'          \code{cns} where there is no data are filled with NA.
#' 
#'          The position indicated by \code{overTox} experience overly toxicity. In the complete single trial, the dose 
#'          level and all the dose levels above will be eliminated.
#'          
#' @return The \code{CFO.next()} function returns a list object comprising the following elements:
#' \itemize{
#'   \item{taget: }{the target DLT rate.}
#'   \item{cys: }{the current counts of DLTs observed at the left, current, and right dose levels}
#'   \item{cns: }{the current counts of patients treated at the left, current, and right dose levels}
#'   \item{decision: }{the decision in the CFO design, where \code{left}, \code{stay}, and \code{right} represent the 
#'   movement directions, and \code{stop} indicates stopping the experiment}
#'   \item{curDoses: }{the current level.}
#'   \item{nextDose: }{the recommended dose level for the next cohort.}
#'   \item{overTox: }{the situation regarding which position experiences overly toxicity, where 'NA' signifies that the 
#'   occurrence of overly toxicity did not happen.}
#' }
#'         
#' @author Jialu Fang and Wenliang Wang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#' 
#' @examples
#' ## determine the dose level for the next cohort of new patients
#' cys <- c(0,1,0); cns <- c(3,6,0)
#' CFO.next(phi=0.2, cys=cys, cns=cns, curDose=3, add.args=list(alp.prior=0.2, bet.prior=0.8))
#' 
#' cys <- c(NA,3,0); cns <- c(NA,3,0)
#' CFO.next(phi=0.2, cys=cys, cns=cns, curDose=1, add.args=list(alp.prior=0.2, bet.prior=0.8))
#' 
#' cys <- c(0,3,NA); cns <- c(3,3,NA)
#' CFO.next(phi=0.2, cys=cys, cns=cns, curDose=7, add.args=list(alp.prior=0.2, bet.prior=0.8))
#' 
#' @import stats
#' @export
CFO.next <- function(phi, cys, cns, curDose,add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, add.args=list()){
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= 0.95) & (add.args$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2)) 
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=1)$value
    const.max <- integrate(fn.max, lower=0, upper=1)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
  }
  
  
  OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
      pC <- 1 - ps$p2
      pL <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL
      
    }else if (type=="R"){
      pC <- 1 - ps$p1
      pR <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }
  
  All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
    for (y1cur in 0:n1){
      for (y2cur in 0:n2){
        ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
      }
    }
    ret.mat
  }
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<phi<upper)
  margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
      dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
  }
  
  # Obtain the table of marginal distribution of (y1, y2) 
  # after intergrate out (phi1, phi2)
  # under H0 and H1
  # H0: phi1=phi, phi < phi2 < 2phi
  # H1: phi2=phi, 0   < phi1 < phi
  margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
      p.y1s <- dbinom(0:n1, n1, phi)
      p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
      p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
      p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
  }
  
  
  # Obtain the optimal gamma for the hypothesis test
  optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H0[1:i])
        err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
    }else if (type=='R'){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H1[1:i])
        err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
      
    }
    list(gamma=gam, min.err=min.err)
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(add.args$alp.prior)){
    add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
  }
  alp.prior <- add.args$alp.prior
  bet.prior <- add.args$bet.prior
  set.seed(seed)
  
  cover.doses <- c(0,0,0)
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.doses[i] <- NA
    }else{
      add.args <- c(list(y=cy, n=cn),list(alp.prior=phi, bet.prior=1-phi))
      if (overdose.fn(phi, add.args)){
        cover.doses[i:3] <- 1
        break()
      }
    }
  }
  cover.doses <- ifelse(is.na(cys), NA, cover.doses)
  
  position <- which(cover.doses == 1)[1]
  overTox <- c("Left", "Current", "Right")[position]
  
  if ((cover.doses[2] == 1)&(curDose == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (cover.doses[2] == 1){
      index <- -1
      decision <- "de-escalation"
    }
    else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        index <- 0
        decision <- "stay"
      }
      else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          index <- -1
          decision <- "de-escalation"
        }else if (!v1 & v2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }
  
  nextDose <- curDose+index
  out <- list(target=phi, cys=cys, cns=cns, decision=decision, curDose = curDose, 
              nextDose=nextDose, overTox=overTox)
  class(out) <- "cfo"
  return(out)
}
