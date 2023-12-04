#' CFO2d.sim.multi
#'
#' This function runs multiple simulations of the 2dCFO design and averages the operating characteristics.
#'
#' @param phi Target dose-intensity level (DIL) rate.
#' @param p.true A matrix representing the true DIL rates under the different dose levels.
#' @param ncohort Integer. The number of cohorts (default is 20).
#' @param cohortsize Integer. The sample size in each cohort (default is 3).
#' @param init.level A numeric vector of length 2 representing the initial dose level (default is c(1,1)).
#' @param add.args An optional list of additional arguments, includes 'alp.prior' (default is phi) and 'bet.prior'(default is 1-phi).
#' @param n_sim Integer. The number of simulations to run (default is 1000).
#' @param seeds A vector of random seed for each simulations, for example, \code{seeds = 1:n_sim} (default is NULL).
#'
#' @return A list with the averaged operating characteristics across all simulations.
#' \itemize{
#'   \item{p.true}{The matrix of the true DLT rates under the different dose levels.}
#'   \item{selPercent}{The matrix of the selection percentage of each dose level}
#'   \item{dose.ns}{A matrix of the averaged number of patients allocated for different doses.}
#'   \item{DLT.ns}{A matrix of the averaged number of DLT observed for different doses.}
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
#' CFO2d.oc(phi=0.3, p.true=p.true, ncohort = 20, cohortsize = 3, n_sim = 1000)

CFO2d.oc <- function(phi, p.true, ncohort = 20, cohortsize = 3, init.level = c(1,1), add.args = list(alp.prior = phi, bet.prior = 1 - phi), seed = NULL, n_sim = 1000) {
  
  # Run the CFO2d.sim function n_sim times using lapply
  results <- lapply(1:n_sim, function(i) {
    CFO2d.sim(phi, p.true, ncohort, cohortsize, init.level, add.args)
  })
  
  selPercent <- matrix(0, dim(p.true)[1], dim(p.true)[2])
  for (i in 1:n_sim) {
    selPercent[results[[i]]$MTD] <- selPercent[results[[i]]$MTD] + 1
  }
  
  # Compute the average of the results
  avg_results <- list()
  avg_results$p.true <- p.true
  avg_results$selPercent <- selPercent / n_sim
  avg_results$dose.ns <- Reduce('+', lapply(results, `[[`, "dose.ns")) / n_sim
  avg_results$DLT.ns <- Reduce('+', lapply(results, `[[`, "DLT.ns")) / n_sim
  avg_results$pcorrect <- mean(sapply(results, `[[`, "correct"))
  avg_results$npercent <- mean(sapply(results, `[[`, "npercent"))
  avg_results$ptoxic <- mean(sapply(results, `[[`, "ptoxic"))
  avg_results$ntox <- mean(sapply(results, `[[`, "ntox"))
  
  return(avg_results)
}


CFO2d.oc <- function(phi, p.true, ncohort = 20, cohortsize = 3, init.level = c(1,1), add.args = list(alp.prior = phi, bet.prior = 1 - phi), n_sim = 1000, seeds = NULL) {
  
  # Run the CFO2d.sim function n_sim times using lapply
  results <- lapply(1:n_sim, function(i) {
    CFO2d.sim(phi, p.true, ncohort, cohortsize, init.level, add.args, seed = seeds[i])
  })
  
  selPercent <- matrix(0, dim(p.true)[1], dim(p.true)[2])
  for (i in 1:n_sim) {
    selPercent[results[[i]]$MTD] <- selPercent[results[[i]]$MTD] + 1
  }
  
  # Compute the average of the results
  avg_results <- list()
  avg_results$p.true <- p.true
  avg_results$selPercent <- selPercent / n_sim
  avg_results$dose.ns <- Reduce('+', lapply(results, `[[`, "dose.ns")) / n_sim
  avg_results$DLT.ns <- Reduce('+', lapply(results, `[[`, "DLT.ns")) / n_sim
  avg_results$pcorrect <- mean(sapply(results, `[[`, "correct"))
  avg_results$npercent <- mean(sapply(results, `[[`, "npercent"))
  avg_results$ptoxic <- mean(sapply(results, `[[`, "ptoxic"))
  avg_results$ntox <- mean(sapply(results, `[[`, "ntox"))
  
  return(avg_results)
}





