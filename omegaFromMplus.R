omegaFromMplus <- function(model, weight=NULL) {
  ## model - Mplus output file
  ## weights - list of vectors with item weights for each group, or a vector with the same weights for all the group, 
  ##           or NULL (all the weights are 1)
  ## WARNING - works only with one factor CFA
  
  require(MplusAutomation)
  model <- readModels(model)
  items <- strsplit(model$input$variable$usevariables, split=" ")[[1]] ## item labels
  n.items <- length(items) ## number of items
  groups <- names(table(model$parameters$unstandardized$Group)) ## group labels
  n.groups <- length(groups)
  parameters <- model$parameters$unstandardized
  parameters$type <- NA
  parameters$type[grep("BY", model$parameters$unstandardized[,1])] <- "BY"
  parameters$type[parameters[,1] == "Intercepts"] <- "INT"
  parameters$type[parameters[,1] == "Residual.Variances"] <- "RES"
  parameters$type[grep("WITH", model$parameters$unstandardized[,1])] <- "WITH"
  parameters$type[parameters[,1] == "Means"] <- "MEAN"
  parameters$type[parameters[,1] == "Variances"] <- "VAR"
  #parameters$first <- NA
  #parameters$first[parameters$type == "WITH"] <- gsub(".WITH", "", parameters["paramHeader", parameters$type == "WITH"])
  omega <- NULL
  
  if (n.groups == 0) {
    groups <- "estimates"
    parameters$Group <- groups
  }
  for (g in groups) {
    if (is.null(weight)) {
      weights.group <- rep(1, n.items)
    } else if (is.vector(weight)) {
      weights.group <- weight[[g]]
    } else {
      weights.group <- weights()
    }
    param <- parameters[parameters$Group == g, ]
    #print(class(param))
    loads <- param[c(1:n.items), 3]
    lat.var <- param[param$type == "VAR", 3]
    res.var <- param[param$type == "RES", 3]
    res.covar <- param[param$type == "WITH", 3]
    
    #print(param)
    #print(g)
    if (n.groups > 0) {
      est.cov <- model$residuals[[g]]$covarianceEst ## estimated covariance matrix
      obs.cov <- model$sampstat[[g]]$covariances ## observed covariance matrix
    } else {
      est.cov <- model$residuals$covarianceEst ## estimated covariance matrix
      obs.cov <- model$sampstat$covariances ## observed covariance matrix
    }
    est.cov[upper.tri(est.cov)] <- t(est.cov)[upper.tri(t(est.cov))] ## fill the upper triangle
    obs.cov[upper.tri(obs.cov)] <- t(obs.cov)[upper.tri(t(obs.cov))] ## fill the upper triangle
    
    
    om1 <- (sum(loads)**2 * lat.var)/(sum(loads)**2 * lat.var + sum(res.var) + sum(res.covar)) ## Raykov (2001) omega
    om2 <- (sum(loads)**2 * lat.var)/sum(weights.group*est.cov) ## Bentler (1972, 2009) omega
    om3 <- (sum(loads)**2 * lat.var)/sum(weights.group*obs.cov) ## McDonald (1999) omega
    alpha <- n.items/(n.items - 1) * (1 - sum(diag(obs.cov))/sum(obs.cov))
    om.group <- c(om1, om2, om3, alpha)
    omega <- cbind(omega, om.group)
  }
  rownames(omega) <- c("Raykov","Bentler","McDonald","alpha")
  colnames(omega) <- groups
  return(omega)
}