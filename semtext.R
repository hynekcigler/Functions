

semtext <- function(fit) {
  fitvalue <- fitMeasures(fit, fit.measures = c("chisq.scaled", "df.scaled", "pvalue.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled", "srmr"))
  if(fitvalue["pvalue.scaled"] < .001) {
    operator1 <- "< .001"
  } else {
    operator1 <- paste0("= ", round(fitvalue["pvalue.scaled"], 3))
  }
  
  result <- paste0("x2(", fitvalue["df.scaled"],") = ", round(fitvalue["chisq.scaled"], 2), ", p ", operator1,
                ", TLI = ", round(fitvalue["tli.scaled"], 3),
                ", RMSEA = ", round(fitvalue["rmsea.scaled"], 3), 
                " (95% CI = [", round(fitvalue["rmsea.ci.lower.scaled"], 3), "; ", round(fitvalue["rmsea.ci.upper.scaled"], 3), "]", 
                # ", p(RMSEA < .05) = ", round(fitvalue["rmsea.pvalue.scaled"], 3), 
                ")", 
                ", SRMR = ", round(fitvalue["srmr"], 3))
  if (!is.na(fitvalue["BIC"])) {
    result <- paste0(result, ", BIC = ", fitvalue["BIC"])
  }
  return(result)
}

semtext2 <- function(fit) {
  fitvalue <- fitMeasures(fit, fit.measures = c("chisq.scaled", "df.scaled", "pvalue.scaled", "tli.robust", "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "rmsea.pvalue.robust", "srmr","bic"))
  if(fitvalue["pvalue.scaled"] < .001) {
    operator1 <- "< .001"
  } else {
    operator1 <- paste0("= ", round(fitvalue["pvalue.scaled"], 3))
  }
  
  result <- paste0("x2(", fitvalue["df.scaled"],") = ", round(fitvalue["chisq.scaled"], 2), ", p ", operator1,
                   ", TLI = ", round(fitvalue["tli.robust"], 3),
                   ", RMSEA = ", round(fitvalue["rmsea.robust"], 3), 
                   " (95% CI = [", round(fitvalue["rmsea.ci.lower.robust"], 3), "; ", round(fitvalue["rmsea.ci.upper.robust"], 3), "]", 
                   # ", p(RMSEA < .05) = ", round(fitvalue["rmsea.pvalue.scaled"], 3), 
                   ")", 
                   ", SRMR = ", round(fitvalue["srmr"], 3))
  if (!is.na(fitvalue["BIC"])) {
    result <- paste0(result, ", BIC = ", fitvalue["BIC"])
  }
  return(result)
}
