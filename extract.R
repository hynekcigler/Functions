semtext <- function(fit, x2=T, cfi=F, tli=T, rmsea=T, rmsea.ci=T, srmr=T, bic=F,
                    clipboard=T, type = "scaled", srmr.type = "srmr", output="text",
                    modelnames=NULL,
                    sep = "\t", dec=",") {
  ## fitted model or list of fitted model
  ## type: standard , scaled , robust
  ## srmr.type: srmr, srmr_bentler, srmr_bentler_nomean, srmr_bollen, srmr_bollen_nomean, srmr_mplus, srmr_mplus_nomean
  
  if(!is.null(modelnames)) {
    model.names <- modelnames
  } else {
    model.names <- as.list(sys.call())[[2]]
    model.names <- gsub("list", "", model.names)[-1]
  }
  
  models <- length(fit)
  fitvalue <- list()
  
  ## extract fitmeasuers
  if (models == 1) {
    fitvalue[[1]] <- fitmeasures(fit)
  } else {
    for (i in c(1:length(fit))) {
      fitvalue[[i]] <- fitmeasures(fit[[i]])
    }
  }

  
  ## type of robust chi2
  if (type == "robust") {
    if(x2 == F) {x2 <- NULL} else {x2 <- c("chisq.scaled", "df.scaled", "pvalue.scaled")}
    if(cfi == F) {cfi <- NULL} else {cfi <- "cfi.robust"}
    if(tli == F) {tli <- NULL} else {tli <- "tli.robust"}
    if(rmsea == F) {rmsea <- NULL} else {rmsea <- "rmsea.robust"}
    if(rmsea.ci == F) {rmsea.ci <- NULL} else {rmsea.ci <- c("rmsea.ci.lower.robust", "rmsea.ci.upper.robust")}
    if(srmr == F) {
      srmr <- NULL
    } else if (is.null(srmr.type)) {
      srmr <- srmr
    } else {srmr <- srmr.type}
    if(bic == F) {bic <- NULL} else {bic <- "bic"}
  }
  else if (type == "scaled") {
    if(x2 == F) {x2 <- NULL} else {x2 <- c("chisq.scaled", "df.scaled", "pvalue.scaled")}
    if(cfi == F) {cfi <- NULL} else {cfi <- "cfi.scaled"}
    if(tli == F) {tli <- NULL} else {tli <- "tli.scaled"}
    if(rmsea == F) {rmsea <- NULL} else {rmsea <- "rmsea.scaled"}
    if(rmsea.ci == F) {rmsea.ci <- NULL} else {rmsea.ci <- c("rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled")}
    if(srmr == F) {
      srmr <- NULL
    } else if (is.null(srmr.type)) {
      srmr <- srmr
    } else {srmr <- srmr.type}
    if(bic == F) {bic <- NULL} else {bic <- "bic"}
  } else {
    if(x2 == F) {x2 <- NULL} else {x2 <- c("chisq", "df", "pvalue")}
    if(cfi == F) {cfi <- NULL} else {cfi <- "cfi"}
    if(tli == F) {tli <- NULL} else {tli <- "tli"}
    if(rmsea == F) {rmsea <- NULL} else {rmsea <- "rmsea"}
    if(rmsea.ci == F) {rmsea.ci <- NULL} else {rmsea.ci <- c("rmsea.ci.lower", "rmsea.ci.upper")}
    if(srmr == F) {
      srmr <- NULL
    } else if (is.null(srmr.type)) {
      srmr <- srmr
    } else {srmr <- srmr.type}
    if(bic == F) {bic <- NULL} else {bic <- "bic"}
  }
  
  ## controls
  if (type == "robust" & is.na(fitvalue[[1]]["cfi.robust"])) {
    stop("Ordered robust statistics do not match estimator. Use different type=")
  } else if (type == "robust" & is.na(fitvalue[[1]]["cfi.robust"])) {
    stop("Ordered scaled statistics do not match estimator. Use different type=")
  }
  
  if (output == "text") {
    result <- NULL
    for (i in c(1:models)) {
      fitvalue1 <- fitvalue[[i]]
      if (!is.null(x2)) {
        if(fitvalue1[x2[3]] < .001) {
          operator1 <- "< 0.001"
        } else {
          operator1 <- paste0("= ", round(fitvalue1[x2[3]], 3))
        }
        x2.out <- paste0("x2(", fitvalue1[x2[2]], ") = ", round(fitvalue1[x2[1]], 1), ", p ", operator1)
      } else {x2.out <- x2}
      if (!is.null(cfi)) {
        cfi.out <- paste0("CFI = ", round(fitvalue1[cfi], 3))
      } else {cfi.out <- cfi}
      if (!is.null(tli)) {
        tli.out <- paste0("TLI = ", round(fitvalue1[tli], 3))
      } else {cfi.out <- tli}
      if (!is.null(rmsea)) {
        rmsea.out <- paste0("RMSEA = ", round(fitvalue1[rmsea], 3))
        if (!is.null(rmsea.ci)) {
          rmsea.out <- paste0(rmsea.out, " with 90% CI = [", 
                              round(fitvalue1[rmsea.ci[1]], 3), ", ", 
                              round(fitvalue1[rmsea.ci[2]], 3), "]")
        }
      } else {rmsea.out <- rmsea}
      if (!is.null(srmr)) {
        srmr.out <- paste0("SRMR = ", round(fitvalue1[srmr], 3))
      } else {srmr.out <- srmr}
      if (!is.null(bic)) {
        bic.out <- paste0("BIC = ", round(fitvalue1[bic], 1))
      } else {bic.out <- bic}
      
      result1 <- paste(c(x2.out, bic.out, cfi.out, tli.out, rmsea.out, srmr.out), collapse = ", ")
      result <- paste0(c(result, result1), collapse = "\n")
    }
    
    
    
    if (clipboard == T) {
      cat(result, file = "clipboard")
      cat("Text output coppied to clipboard: \n\n")
    }
    return(cat(result))
    
    
    
  } else if (output == "table") {
    if (!is.null(x2)) {
      df <- sapply(fitvalue, '[[', x2[2])
      p <- round(sapply(fitvalue, '[[', x2[3]), 3)
      x2 <- round(sapply(fitvalue, '[[', x2[1]), 1)
      x2 <- cbind(x2, df, p)
    }
    
    if (!is.null(cfi)) {
      CFI <- round(sapply(fitvalue, '[[', cfi), 3)
    } else {CFI <- cfi}
    
    if (!is.null(tli)) {
      TLI <- round(sapply(fitvalue, '[[', tli), 3)
    } else {TLI <- tli}
    
    if (!is.null(rmsea)) {
      RMSEA <- round(sapply(fitvalue, '[[', rmsea), 3)
      if (!is.null(rmsea.ci)) {
        lower <- round(sapply(fitvalue, '[[', rmsea.ci[1]), 3)
        upper <- round(sapply(fitvalue, '[[', rmsea.ci[2]), 3)
      } else {lower <- upper <- rmsea.ci}
      RMSEA <- cbind(RMSEA, lower, upper)
    } else {RMSEA <- rmsea}
    
    if (!is.null(tli)) {
      SRMR <- round(sapply(fitvalue, '[[', srmr), 3)
    } else {SRMR <- srmr}
    
    if (!is.null(bic)) {
      BIC <- round(sapply(fitvalue, '[[', bic), 1)
    } else {BIC <- bic}
    
    result <- cbind(x2, BIC, CFI, TLI, RMSEA, SRMR)
    rownames(result) <- model.names
    if(clipboard == T) {
      cat("Output coppied to clipboard: \n\n")
      write.table(result, "clipboard", sep = sep, dec = dec, col.names = NA, quote = T)
    }
    return(result)
  }
}


semtext2 <- function(fit, x2=T, cfi=F, tli=T, rmsea=T, rmsea.ci=T, srmr=T, bic=F,
                     clipboard=T, type = "robust", srmr.type = "srmr", output="text",
                     modelnames=NULL,
                     sep = "\t", dec=",") {
  ## only abbreviation function for type=robust
  return(semtext(fit=fit, x2=x2, cfi=cfi, tli=tli, rmsea=rmsea, rmsea.ci=rmsea.ci, srmr=srmr, bic=bic,
                 clipboard=clipboard, type = type, srmr.type = srmr.type, output=output,
                 modelnames=modelnames,
                 sep = sep, dec=dec))
}
