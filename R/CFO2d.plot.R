

plot_dose <- function(CFO2d.res){
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


# Usage

CFO2d.res <- CFO2d.sim(phi=0.3, p.true, ncohort = 20, cohortsize = 3)
plot_dose(CFO2d.res)
