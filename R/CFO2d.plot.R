#' Plotting Dose Function
#'
#' This function takes the result of CFO2d.sim simulation and plots the combined dose level 
#' for each sequence of patients treated. The doses are displayed as a step plot, 
#' with points representing each patient. The point is filled if a dose-limiting toxicity (DLT) is observed.
#'
#' @param CFO2d.res A list result from the CFO2d.sim function. It should include the following components:
#' \itemize{
#'  \item{dose}{A matrix of doses given to each cohort.}
#'  \item{DLT}{A vector of DLTs observed in each cohort.}
#'  \item{dose.ns}{A vector of the number of subjects at each dose.}
#'  \item{p.true}{A matrix representing the true probabilities of DLT at each dose combination.}
#' }
#'
#' @return A ggplot object of the dose assigned to each patient, and whether DLT was observed.
#' The x-axis represents the sequence of patients treated, while the y-axis represents the combined dose level.
#' The points are filled if a DLT is observed.
#' @import ggplot2
#' @export
#' @examples
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#' CFO2d.res <- CFO2d.sim(phi=0.3, p.true, ncohort = 20, cohortsize = 3)
#' CFO2d.plot(CFO2d.res)


CFO2d.plot <- function(CFO2d.res){
  dose <- CFO2d.res$dose
  DLT <- CFO2d.res$DLT
  ncohort <- length(CFO2d.res$DLT)
  cohortsize <- sum(CFO2d.res$dose.ns)/ncohort
  dim <- dim(CFO2d.res$p.true)
  
  # Generate y_labels
  y_labels <- expand.grid(1:dim[1], 1:dim[2])
  y_labels <- apply(y_labels, 1, function(x) paste('(', x[1], ',', x[2], ')'))
  
  # Generate sequences for each patient
  sequences <- 1:(ncohort * cohortsize)
  
  # Generate dose_levels for each patient
  dose_levels <- rep(match(apply(dose, 1, function(x) paste('(', x[1], ',', x[2], ')')), y_labels), each = cohortsize)
  
  # Generate DLT_observed for each patient
  DLT_observed <- unlist(mapply(function(dlt, size) c(rep(1, dlt), rep(0, size - dlt)), DLT, rep(cohortsize, ncohort)))
  
  df <- data.frame(sequence = sequences, dose_level = dose_levels, DLT_observed = DLT_observed)
  
  # Create the plot
  p <- ggplot(df, aes(x = sequence, y = dose_level)) +
    geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
    geom_step(direction = 'hv', color = 'black') +
    scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
    labs(x = "Sequence of patients treated", 
         y = "Combined dose level",
         fill = 'DLT observed') +
    theme_minimal() +
    theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
    scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
  
  # Display the plot
  print(p)
}

