closefit <- function(m0, m1, df.fix=0) {
  ## This function compute indexes of the close fit of IRT models produced by the lme4 packages.
  ## m0, m1 are the the object of the class merMod
  ## df.fix = glmer function needs to have at least one random effect. If there is no random effect
  ##          in null model, you have to include "phantom" random effect alongside the fixed effect,
  ##          e.g. resp ~ 1 + items + (1 | items) . In such, case, set df.fix =1
  ## NOTE: The indices are not based on M2 statistics as was suggested 
  ##       by Maydeu-Olivares and Joe (2006) or Cai and Hansen (2013), 
  ##       but it is based purely on the model deviance. This indices
  ##       are not comparable to the traditional fit indices in CFA or 
  ##       in confirmatory factor models.
  ## m0 = baseline model
  ## m1 = research model
  
  cat("NOTE: The indices are not based on M2 statistics as was suggested by Maydeu-Olivares and Joe (2006) or Cai and Hansen (2013), but it is based purely on the model deviance. This indices are not comparable to the traditional fit indices in CFA or in confirmatory factor models.
      
")
  
  m0 <- summary(m0)
  m1 <- summary(m1)
  cfi <- ((m0$AICtab["deviance"]-(m0$AICtab["df.resid"]-df.fix))-
            (m1$AICtab["deviance"]-m1$AICtab["df.resid"]))/
    (m0$AICtab["deviance"]-(m0$AICtab["df.resid"]-df.fix))
  tli <- ((m0$AICtab["deviance"]/(m0$AICtab["df.resid"]-df.fix))-
            (m1$AICtab["deviance"]/m1$AICtab["df.resid"]))/
    (m0$AICtab["deviance"]/(m0$AICtab["df.resid"]-df.fix)-1)
  if (m1$AICtab["deviance"] - m1$AICtab["df.resid"] < 0) {
    rmsea <- NaN
    warning("Chi2 is less than number of df - perfect fit. Fit indices could be biased", call.=F)
  } else {
    rmsea <- sqrt(m1$AICtab["deviance"] - m1$AICtab["df.resid"])/
      sqrt(m1$AICtab["df.resid"] * (m1$devcomp$dims["N"] - 1))
  }
  
  cfi <- max(min(cfi, 1), 0)
  cfi <- round(as.numeric(cfi), 3)
  tli <- round(as.numeric(tli), 3)
  rmsea <- round(as.numeric(rmsea), 3)
  
  chi2 <- round(as.numeric(m1$AICtab["deviance"]), 3)
  df <- round(as.numeric(m1$AICtab["df.resid"]), 3)
  p <- round(as.numeric(pchisq(chi2, df, lower.tail = F)), 3)
  
  result <- c(chi2=chi2, df=df, p=p, cfi=cfi, tli=tli, rmsea=rmsea)

  return(result)
}
