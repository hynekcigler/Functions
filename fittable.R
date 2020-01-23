fittable <- function(..., indices="robust", all=F) {
  require(lavaan)
  models <- list(...)
  indic.names <- c("x2", "df", "p", "BIC", "CFI", "TLI", "RMSEA", "lower.RMSEA", "upper.RMSEA", "SRMR")
  
  if (indices == "scaled") {
    indic <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "bic", "cfi.scaled", "tli.scaled", "rmsea.scaled",
               "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr")

  } else if (indices == "robust") {
    indic <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "bic", "cfi.robust", "tli.robust", "rmsea.robust",
               "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr")
  } else {
    indic.names <- indic <- indices
  }

  table.fit <- NULL
  all.fit <- list()
  for (m in c(1:length(models))) {
    all.fit[[m]] <- fitMeasures(models[[m]])
    fits2 <- all.fit[[m]][indic]
    table.fit <- rbind(table.fit, fits2)
  }
  rownames(table.fit) <- paste("model", c(1:length(models)))
  colnames(table.fit) <- indic.names
  names(all.fit) <- paste0("model", c(1:length(models)))
  
  print(round(table.fit, digits=3))
  if (all == F) {
    return(invisible(table.fit))
  } else {
    return(invisible(list(table=table.fit, all=all.fit)))
  }
  
  
}
