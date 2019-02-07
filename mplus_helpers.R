#****** Programmer 1  ©:        Hynek Cígler                                           ******# 
#****** Programmer 2  ©:                                                               ******#
#****** Contact/Email:          hynek.cigler@mail.muni.cz                              ******#
#****** Affiliation:            Masaryk University                                     ******#
#****** Project:                TALIS MS 2018                                          ******#
#****** Analysis:               Scaling TALIS MS                                       ******#
#****** Last Update:            September 2018                                         ******#
##############################################################################################




# Reliability of factor score estimates -----------------------------------

## This functions compute reliability of factor score estimates using the empirical data.
## model   = object produced by function readModels or the path to the Mplus output (only one file/Mplus object supported)
## scale   = name of scale or vector of scale names in the case of multidimensional CFA (the same as in the Mplus input)
## country = if result should be provided for each country separately, name of country variable
## med     = provide also median across countries?

empirical_rxx <- function(model, scale=NULL, country=NULL, med=F) {
  model <- model.reading(model, recursive = F, what = c("parameters", "savedata"))
  if (is.null(model$savedata)) {
    stop("Model does not contain data.")
  } else if (model$input$analysis$estimator == "WLSMV")  {
    warning("Estimation is not possible with categorical model estimated by WLSMV estimator.")
    return(NA)
  } else {
    dat <- model$savedata ## data with factor scores
  }
  
  if(is.null(scale)) {
    scale <- names(table(model$parameters$unstandardized[model$parameters$unstandardized$paramHeader == "Variances" , "param"]))
  }
  
  err <- paste0(scale, "_SE")      ## name of error variables
  
  rxx <- NULL
  
  if(is.null(country)) {
    for (i in c(1:length(scale))) {
      rmse <- sqrt(mean(dat[[err[i]]]**2)) ## RMSE
      truevar <- var(dat[[scale[i]]])      ## variance of factor scores
      rxx1 <- truevar/(truevar + rmse**2)  ## scale reliability
      rxx <- c(rxx, rxx1)                  ## merge results
    }
    names(rxx) <- scale                    ## name results
  } else {
    
    countryvar <- country
    country <- names(table(model$savedata[[countryvar]]))
    
    for (g in c(1:length(country))) {
      rxx2 <- NULL
      for (i in c(1:length(scale))) {
        rmse <- sqrt(mean(dat[dat[[countryvar]] == country[g], err[i]]**2)) ## RMSE
        truevar <- var(dat[dat[[countryvar]] == country[g], scale[i]])      ## variance of factor scores
        rxx1 <- truevar/(truevar + rmse**2)                                 ## scale reliability
        rxx2 <- c(rxx2, rxx1)                                               ## merge results
      }
      rxx <- rbind(rxx, rxx2)                                               ## merge results
    }
    colnames(rxx) <- scale; rownames(rxx) <- country
    
    if (med == T) {
      rxx <- list(country = rxx, median = apply(rxx, 2, median))
    }
  }
  return(rxx)
}



# Identification of Heywood case ------------------------------------------

heywood <- function(model) {
  model <- model.reading(model, recursive = F, what = c("parameters"))
  negative <-   sum(model$parameters$unstandardized$CONFIGURAL.MODEL$est[
    model$parameters$unstandardized$CONFIGURAL.MODEL$paramHeader == "Residual.Variances"
    ] <= 0)
  if (negative > 0) {
    warning("Heywood case!")
  }
  return(negative)
}



# Modification indices extraction -----------------------------------------

## model = path to folder, Mplus output, or result of readModels function
## n = if null, all M.I.s showed; otherwise n-first M.I.s showed (must be sort.=T)
## cut = show only M.I.s with chi2 > cut
## sort. = should be output sorted by M.I.s?
## sort.by = if sort.=T, then what is the sorter. Applicable headings of the result - "MI", "EPC", "Std_EPC" or "StdYX_EPC".
## recursive = if model is the path to a folder, should the function look also into subfolders?

## Return (sorted, subseted) data.frame with M.I.s, or NULL and warning if there are not any M.I.s to display
## In the case of more models, return list of above mentioned results

modindices <- function(model, n = NULL, cut=3.84, sort. = T, sort.by = "MI", recursive = F) {
  model <- model.reading(model, recursive = recursive, what="mod_indices")
  if(class(model)[[1]] == "mplus.model") {
    return(modindices1(model, n = n, cut = cut, sort. = sort., sort.by = sort.by))
  } else {
    models <- length(model)
    result <- list()
    try(for (i in c(1:models)) {
      result[i] <- list(modindices1(model[[i]], n = n, cut = cut, sort. = sort., sort.by = sort.by))
      names(result)[i] <- model[[i]]$summaries$Title
    })
    return(result)
  }
}




# Model fit ---------------------------------------------------------------

## Only wraper for SummaryTable function from MplusAutomation package. 
## Function has preseted fit indices and it shorts the column labels

## model = path to folder, Mplus output, or result of readModels function
## info = should be Mplus version, estimator, number of observed variables and groups and sample size printed too?
## recursive = if model is the path to a folder, should the function look also into subfolders?
## fits = prespecified model fits
## nullrmsea = should be also RMSEA for baseline model printed?

modelfit <- function(model, recursive=F, fits = c("Title", 
                                                  "ChiSqM_Value", "ChiSqM_DF", "ChiSqM_PValue", "CFI", "TLI", 
                                                  "RMSEA_Estimate", "RMSEA_90CI_LB", "RMSEA_90CI_UB", "SRMR", "WRMR"),
                     info = T, nullrmsea=F, clipboard = F, dec=",") {
  require(MplusAutomation)
  require(ddpcr)
  if (info == T) {
    fits <- c(fits, "Mplus.version", "Estimator", "Observations", "NGroups", "NDependentVars")
  }
  model <- model.reading(model, recursive = recursive, what = c("input", "summaries"))

  quiet(result <- SummaryTable(model, include.rownames = T, keepCols = fits, display = F))
  names(result)[names(result) == "ChiSqM_Value"] <- "x2"
  names(result)[names(result) == "ChiSqM_DF"] <- "df"
  names(result)[names(result) == "ChiSqM_PValue"] <- "p"
  names(result)[names(result) == "RMSEA_Estimate"] <- "RMSEA"
  names(result)[names(result) == "RMSEA_90CI_LB"] <- "CI_l"
  names(result)[names(result) == "RMSEA_90CI_UB"] <- "CI_u"
  names(result)[names(result) == "Observations"] <- "N"
  names(result)[names(result) == "NDependentVars"] <- "items"
  names(result)[names(result) == "NGroups"] <- "groups"
  names(result)[names(result) == "Mplus.version"] <- "version"
  if (nullrmsea == T) {
    result$nullRMSEA <- round(nullRMSEA(model)[,2], 3)
  }
  if (clipboard == T) {
    cat("Results were copied to clipboard.")
    write.table(result, "clipboard", col.names = F, sep="\t", dec=dec)
  }
  
  return(result)
}


# RMSEA for the null model ------------------------------------------------

## model = path to folder, Mplus output, or result of readModels function
## recursive = if model is the path to folder, should the function look also to the subfolders?

## In the case of one model, return RMSEA of its null model.
## In the case of more models, return matrix of model titles and their null models.

nullRMSEA <- function(model, recursive=F) {
  model <- model.reading(model, recursive = recursive) ## just read model/models
  
  if(class(model)[[1]] == "mplus.model") {
    nullRMSEA1(model)
  } else {
    models <- length(model)
    result <- data.frame()
    for (i in c(1:models)) {
      result[i, 1] <- model[[i]]$summaries$Title
      result[i, 2] <- nullRMSEA1(model[[i]])
    }
    names(result) <- c("model", "rmsea")
    result$low <- result$rmsea < .158
    return(result)
  }
}


######################################################################################################

# Support function for the above functions -----------------------------------------------------------
# Should not be used directly


## Only read models
model.reading <- function(model, recursive=F, what="all") {
  require(MplusAutomation) ## load package if needed
  
  if ((class(model)[1] != "mplus.model") & (class(model)[1] != "mplus.model.list")) { ## read file if model is path to file
    model <- readModels(model, recursive=recursive, what=what)
  }
  return(model)
}

## Just compute RMSEA of the null model
nullRMSEA1 <- function(model) {
  x2 <- model$summaries$ChiSqBaseline_Value ## null chi2
  df <- model$summaries$ChiSqBaseline_DF ## null df
  G <- model$summaries$NGroups ## number of groups
  N <- model$summaries$Observations ## sample size
  if (is.null(N)) {
    stop("Multigroup analysis is not supported.")
  }
  
  rmsea <- sqrt(x2/((N-1)*df) - 1/(N-1))*sqrt(G) ## compute null RMSEA
  
  if (rmsea < .158) {
    warning("RMSEA of the null model is lower than 0.158 suggested by Kenny, Kaniskan, & McCoach (2015). 
            Incremental fit indices such as TLI, CFI would not be informative.")
  }
  return(rmsea)
}

## Just exctact fit indices
modindices1 <- function(model, n = NULL, cut=3.84, sort. = T, sort.by = "MI") {
  if(is.null(model$mod_indices)) {
    warning("No modification indices available for the model.")
    return(NULL)
  }
  MI <- model$mod_indices
  MI <- MI[MI$MI > cut, ]
  if(nrow(MI) == 0) {
    warning("No modification indices above the cut value.")
    return(NULL)
  }
  if(nrow(MI) == 1 & MI$MI[1] == 999) {
    warning("Modification indices are not meaningfull - model has perfect fit (probably only three items).")
    return(NULL)
  }
  if (sort. == T) {
    MI <- MI[order(abs(MI[[sort.by]]), decreasing = T), ]
  }
  if (!is.null(n) & sort. == T) {
    return(head(MI, n = n))
  } else {
    return(MI)
  }
}
