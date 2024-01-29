#' This function runs multiple simulations of the 2dCFO design and averages the operating characteristics.
#'
#' @param target the target DLT rate.
#' @param p.true a matrix representing the true DIL rates under the different dose levels.
#' @param ncohort the total number of cohorts, the default value is 20.
#' @param cohortsize the number of patients or size of each cohort. the default value is 3. 
#' @param init.level a numeric vector of length 2 representing the initial dose level (default is c(1,1)).
#' @param prior.para the prior parameters for a beta distribution, usually set as list(alp.prior=target, bet.prior=1-target) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param nsimu the total number of trials to be simulated. The default value is 1000.
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli = 0.95}) for general use.
#' @param extrasafe set \code{extrasafe = TRUE} to impose a more strict early stopping rule for
#'                   extra safety.
#' @param offset a small positive number (between \code{0} and \code{0.5}) to control how strict the
#'                stopping rule is when \code{extrasafe=TRUE}. A larger value leads to
#'                a more strict stopping rule. The default value \code{offset = 0.05}
#'                generally works well.
#' @param seeds A vector of random seed for each simulations, for example, \code{seeds = 1:nsimu} (default is NULL).
#'
#' @return A list with the averaged operating characteristics across all simulations.
#' \itemize{
#'   \item{p.true}{The matrix of the true DLT rates under the different dose levels.}
#'   \item{selPercent}{The matrix of the selection percentage of each dose level}
#'   \item{dose.ns}{A matrix of the averaged number of patients allocated for different doses in one simulation.}
#'   \item{DLT.ns}{A matrix of the averaged number of DLT observed for different doses in one simulation.}
#'   \item{pcorrect}{The percentage of the correct selection of the MTD.}
#'   \item{npercent}{The averaged percentage of patients assigned to the target DLT rate.}
#'   \item{ptoxic}{The averaged percentage of patients assigned to dose levels with a DLT rate greater than the target.}
#'   \item{ntox}{The averaged total number of DLTs observed.}
#' }
#' @export
#' @examples
#' ## Simulate a two-dimensional dose-finding trial with 20 cohorts of size 3 for 1000 replications.
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#' 0.10, 0.15, 0.30, 0.45, 0.55,
#' 0.15, 0.30, 0.45, 0.50, 0.60), 
#' nrow = 3, ncol = 5, byrow = TRUE)
#' 
#' CFO2doc <- CFO2d.oc(target=0.3, p.true=p.true, ncohort = 25, cohortsize = 3, seeds = 1:10, nsimu = 10)
#' summary(CFO2doc)
#' plot(CFO2doc)

CFO2d.oc <- function(target, p.true, ncohort = 20, cohortsize = 3, init.level = c(1,1), prior.para = list(alp.prior = target, bet.prior = 1 - target), cutoff.eli=0.95, extrasafe=FALSE, offset=0.05, nsimu = 1000, seeds = NULL) {
  
  # Run the CFO2d.simu function nsimu times using lapply
  results <- lapply(1:nsimu, function(i) {
    CFO2d.simu(target, p.true, ncohort, cohortsize, init.level, prior.para, seed = seeds[i], cutoff.eli=cutoff.eli, extrasafe=extrasafe, offset=offset)
  })
  
  selPercent <- matrix(0, dim(p.true)[1], dim(p.true)[2])
  for (i in 1:nsimu) {
    selPercent[results[[i]]$MTD] <- selPercent[results[[i]]$MTD] + 1
  }
  
  # Compute the average of the results
  avg_results <- list()
  avg_results$p.true <- p.true
  avg_results$selPercent <- selPercent / nsimu
  avg_results$dose.ns <- Reduce('+', lapply(results, `[[`, "dose.ns")) / nsimu
  avg_results$DLT.ns <- Reduce('+', lapply(results, `[[`, "DLT.ns")) / nsimu
  avg_results$pcorrect <- mean(sapply(results, `[[`, "correct"))
  avg_results$npercent <- mean(sapply(results, `[[`, "npercent"))
  avg_results$ptoxic <- mean(sapply(results, `[[`, "ptoxic"))
  avg_results$ntox <- mean(sapply(results, `[[`, "ntox"))
  
  class(avg_results) <- "cfo"
  
  return(avg_results)
}
