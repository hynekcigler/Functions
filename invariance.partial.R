invariance.partial <- function(model, type, MI.cut=5, min.anchor=3, write.model=T, write.pattern=T, name.model=NULL, 
                               epc.cut=0, order.by="mi", .order=F) {
  ## model = Mplus output of metric or scalar model; alternatively object of class mplus.model (result of function readModels)
  ## type = c("metric", "scalar")
  ## min.anchor - can be number of minimum number of anchored items for each group, or it can be "all" which select or proposed items with M.I. < MI.cut
  ## order.by - drop-out items by M.I. (order.by="MI") or E.P.C. (order.by="epc")? 
  ## Mi.cut - minimum M.I. value of items which are suggested to be freely estimate
  ## epc.cut - minimum value of xy standardized E.P.C. change (expected parameter change under xystd)
  ## write.model - should. be the updated part of model printed?
  ## name.model - name of the file into which the updated model are saved. If null, model will be printed only in consol.
  ## .order - if order = T, item M.I. output is printed in the M.I. (if order.by="mi") or E.P.C. order (if order.by="epc"). 
  ##          If order=F, item are printed in the same order as they are listed in usevariable parameter in Mplus (by group)
  ## write.pattern - if TRUE, pattern of suggested items to set free are printed together with tests of dependences (Fleiss kappa). 
  ##                 (Do not trust these tests.)
  
  require(MplusAutomation)
  
  if ((class(model) != "mplus.model")[1]) {
    model <- readModels(model)
  }
  
  items <- strsplit(model$input$variable$usevariables, split=" ")[[1]] ## item labels
  n.items <- length(items) ## number of items
  groups <- names(table(model$parameters$unstandardized$Group)) ## group labels
  
  if (min.anchor == "all") { ## if "all", all items with MI and std E.P.C higher than cuts are printed
    min.anchor <- 0
  } else if (n.items < min.anchor) { ## if number of anchored items is to high, it uses number of items minus 2
    warning("Minimum value of anchoring items (min.anchor = ", min.anchor, ") is lower than the actual number of items at least in one group. 
            Minimum number of anchored items was set to ", n.items-2, ".", call. = F)
    min.anchor <- n.items-2
  }
  max.items <- n.items - min.anchor ## maximum number of items which can be released
  
  MI <- model$mod_indices ## extract modification indices
  MI2 <- MI 
  
  text.model <- "" ## just prepare empty object for subsequent use
  
  if (type == "metric") {
    MI2 <- MI2[MI2$operator == "BY" & MI2$MI > MI.cut  & (abs(MI2$StdYX_EPC) > epc.cut),] ## M.I. which meet criteria (M.I. and std. E.P.C. cutt-offs)
    MI3 <- NULL
    for (g in groups) {
      MIx <- MI2[(MI2$Group == g), ]
      if (order.by == "mi") {
        MIx <- MIx[(MIx$MI %in% tail(MIx$MI, max.items)), ] ## select items by MI
      } else if (order.by == "epc") {
        MIx <- MIx[(MIx$StdYX_EPC %in% tail(MIx$StdYX_EPC, max.items)), ] ## select items by EPC
      } else {
        stop("Unknown 'order.by' parameter.")
      }
      if (nrow(MIx) == 0) next ## skip if there are now suggested items to release for particular group
      MI3 <- rbind(MI3, MIx)
      if (write.model == T) {
        model.group <- paste("MODEL", g,": ") ## group model header
        for (i in c(1:nrow(MIx))) {
          model.part <- paste(paste(c(MI3[i,1:3]), collapse=" ", sep=""), "*; ", collapse = "", sep="") ## submodel for particular item
          model.group <- paste(model.group, model.part, sep="\n") ## merge with previous items
        }
        text.model <- paste(text.model, model.group, sep="\n\n") ## merge with other groups
      }
    }
  } else if (type == "scalar") {
    MI2 <- MI2[MI2$operator == "NA" & MI2$MI > MI.cut & (abs(MI2$StdYX_EPC) > epc.cut),] ## M.I. which meet criteria (M.I. and std. E.P.C. cutt-offs)
    MI3 <- NULL
    MI3 <- NULL
    for (g in groups) {
      MIx <- MI2[(MI2$Group == g), ]
      if (order.by == "mi") {
        MIx <- MIx[(MIx$MI %in% tail(MIx$MI, max.items)), ] ## select items by MI
      } else if (order.by == "epc") {
        MIx <- MIx[(MIx$StdYX_EPC %in% tail(MIx$StdYX_EPC, max.items)), ] ## select items by EPC
      } else {
        stop("Unknown 'order.by' parameter.")
      }
      if (nrow(MIx) == 0) next ## skip if there are now suggested items to release for particular group
      MI3 <- rbind(MI3, MIx)
      if (write.model == T) {
        model.group <- paste("MODEL", g,": ") ## group model header
        for (i in c(1:nrow(MIx))) {
          model.part <- paste(c(MI3[i,1], "; "), collapse=" ")
          model.group <- paste(model.group, model.part, sep="\n") ## merge with previous items
        }
        text.model <- paste(text.model, model.group, sep="\n\n") ## merge with other groups
      }
    }
    text.model <- gsub("]", "*]", text.model) ## "release" parameters
  } else {
    stop("Error: Unknown type of partial invariance model.") ## e.g. strict model are not supported yet
  }
  
  if (is.null(MI3)) {
    stop("There is now suggested parameters to release. Function has to stop.")
  }
  
  if (write.model == T) {
    if (!is.null(name.model)) {
      cat(text.model, file=name.model) ## write file with text, which can be copy-past to Mplus input
    }
  } else {
    text.model <- NULL
  }
  
  if (write.pattern == F) {
    text.pattern <- NULL
    p.pattern.item <- NULL
    p.pattern.group <- NULL
  } else {
    if (type == "scalar") { ## prepare item column for scalar models
      MI3$item <- MI3$modV1
      MI3$item <- gsub("]", "", MI3$item)
      MI3$item <- gsub("\\[", "", MI3$item)
    } else { ## prepare item column for metric models
      MI3$item <- MI3$modV2
    }
    text.pattern <- matrix(NA, ncol=n.items, nrow=length(groups), dimnames = list(groups, items)) ## prepare pattern matrix
    for (i in c(1:length(groups))) {
      for (j in c(1:n.items)) {
        text.pattern[i,j] <- as.numeric(items[j] %in% MI3[MI3$Group == groups[i],"item"]) ## impute pattern matrix. 1 = shoudl be free, 0 = should be constrained
      }
    }
    text.pattern <- as.data.frame(text.pattern) ## just for subsequent analyses
    if (require(irr, quietly = T)) {
      p.pattern.group <- kappam.fleiss(t(text.pattern)) ## are the items treated by groups similarly? Non-significant value -> OK (no relations)
      p.pattern.group <- c(kappa=p.pattern.group$value, z=p.pattern.group$statistic, p.value=p.pattern.group$p.value)
      p.pattern.item <- kappam.fleiss(text.pattern) ## are the groups treated by items similarly? Non-significant value -> OK (no relations)
      p.pattern.item <- c(kappa=p.pattern.item$value, z=p.pattern.item$statistic, p.value=p.pattern.item$p.value)
    } else {
      p.pattern.group <- p.pattern.item <- "irr package is missing. Cannot compute probability of random pattern in anchor items."
      warning(p.pattern.group)
    }
  }
  if (.order == T) {
    if (order.by == "mi") {
      MI3 <- MI3[order(-MI3$MI),]
    } else {
      MI3 <- MI3[order(-abs(MI3$StdYX_EPC)),]
    }
  }
  
  result <- list(release=MI3, model=cat(paste(text.model, "\n\n")), pattern=text.pattern, p.pattern.group=p.pattern.group, 
                 p.pattern.item=p.pattern.item)
  return(result)
}