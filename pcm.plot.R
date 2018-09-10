### (c) Hynek Cígler | Masaryk University | 2018           ###
### hynek.cigler@mail.muni.cz | http://psych.fss.muni.cz   ###


pcm.plot <- function(theta=seq(-5,5, by=.001), b=0, thresholds=c(-1,1), show="none", ylim=NULL, 
                     color = T, from=0, to=.6) {
  ## theta = just a vector for subsequent ability values, for which the diagram should be plotted
  ## b = item difficulty
  ## thresholds = vector of thresholds. Its sum has to be 0.
  ## show = type of plot
  ##        trace - traceline (item cattegory characteristic curves)
  ##        score - scoring function (expected raw score)
  ##        info - information function (both cattegories and the whole item)
  ##        color - if TRUE, plot colors; if FALSE, use gray scale
  ##        from, to - parameters if color=TRUE
  
  
  if(abs(sum(thresholds)) >.0000001) {
    stop("Sum of thresholds is not zero!")
  }
  
  
  # Category probability scores ---------------------------------------------
  
  pred <- matrix(nrow = length(theta), ncol = length(thresholds)+1, NA)
  
  ## jmenovatel
  x3 <- rep(0, length(theta))
  for (j in c(1:length(thresholds))) {
    x2 <- rep(0, length(theta))
    for (k in c(1:j)) {
      x <- theta-(b-thresholds[k])
      x2 <- x2+x
    }
    x3 <- x3+exp(x2)
  }
  
  pred[,1] <- 1/(1 + x3) ## kat 0
  
  ## citatel
  for (k in c(1:length(thresholds))) {
    x4 <- rep(0, length(theta))
    for(obs in c(1:k)) {
      x <- theta-(b-thresholds[obs])
      x4 <- x+x4
    }
    pred[,1+k] <- exp(x4)/(1+x3)
  }
  
  xij <- t(t(pred)*(c(1:ncol(pred))-1))
  score <- rowSums(xij)
  a <- c(1:(ncol(pred)-1))
  rasch_sc <- rasch <- half <- rep(NA, length(a))
  
  
  # Tresholds ---------------------------------------------------------------
  
  for (i in a) {
    half[i] <- theta[which(abs(score-(a[i]-.5))==min(abs(score-(a[i]-.5))))]
    if (i == 1) {
      rasch[i] <- theta[which(abs(pred[, 1] - rowSums(pred[,c(i:max(a))+1])) == min(abs(pred[, 1] - rowSums(pred[,c(i:max(a))+1]))))]
      rasch_sc[i] <- score[which(abs(pred[, 1] - rowSums(pred[,c(i:max(a))+1])) == min(abs(pred[, 1] - rowSums(pred[,c(i:max(a))+1]))))]
    } else if (i == max(a)) {
      rasch[i] <- theta[which(abs(rowSums(pred[,c(1:i)]) - pred[,max(a)+1]) == min(abs(rowSums(pred[,c(1:i)]) - pred[,max(a)+1])))]
      rasch_sc[i] <- score[which(abs(rowSums(pred[,c(1:i)]) - pred[,max(a)+1]) == min(abs(rowSums(pred[,c(1:i)]) - pred[,max(a)+1])))]
    } else {
      rasch[i] <- theta[which(abs(rowSums(pred[,c(1:i)]) - rowSums(pred[,c(i:max(a))+1])) == min(abs(rowSums(pred[,c(1:i)]) - rowSums(pred[,c(i:max(a))+1]))))]
      rasch_sc[i] <- score[which(abs(rowSums(pred[,c(1:i)]) - rowSums(pred[,c(i:max(a))+1])) == min(abs(rowSums(pred[,c(1:i)]) - rowSums(pred[,c(i:max(a))+1]))))]
    }
    # 
  }
  
  
  # Information -------------------------------------------------------------
  
  ## all
  x <- pred
  for (i in c(1:ncol(pred))) {
    x[,i] <- ((i-(score+1))**2 * pred[,i])
  }
  suminfo <- rowSums(x)
  
  
  ## category
  info <- pred*matrix(as.vector(suminfo), nrow=length(suminfo), ncol = ncol(pred))
  
  # Plots --------------------------------------------------------------------
  
  if (isTRUE(color)) {
    cols <- rainbow(ncol(pred))
  } else {
    cols <- gray(seq(from, to, length.out = ncol(pred)))
  }
  
  
  if(show == "trace") {
    plot(theta, pred[,1], type = "l", lwd=3, col=cols[1], xlab="trait", ylab = "probability", main="Category probability functions")
    for (i in c(2:ncol(pred))) {
      lines(theta, pred[,i], type = "l", lwd=3, col=cols[i])
    }
    if(isTRUE(color)) {
      abline(v=half, col="green", lty=2)
      abline(v=rasch, col="blue", lty=2)
    } else {
      abline(v=half, col="black", lty=3)
      abline(v=rasch, col="black", lty=4)
    }
    if(isTRUE(color)) {
      legend("right", legend=c(c(1:ncol(pred))-1, "half-point", "Rasch-Andrich"), col=c(cols, "green", "blue"), 
             lwd=c(rep(3, ncol(pred)), 1, 1), lty=c(rep(1, ncol(pred)), 2, 2), inset=.01, bg = "white")
      
    } else {
      legend("right", legend=c(c(1:ncol(pred))-1, "half-point", "Rasch-Andrich"), col=c(cols, "black", "black"), 
             lwd=c(rep(3, ncol(pred)), 1, 1), lty=c(rep(1, ncol(pred)), 3, 4), inset=.01, bg = "white")
      
    }
  } else if (show == "score") {
    plot(theta, score, type = "l", lwd=3, col="black", xlab="trait", ylab = "expected score", main="Expected score function")
    # abline(v=half, lty=2, col="gray")
    for (i in c(1:length(half))) {
      if(isTRUE(color)) {
        # half-tresholds
        lines(c(half[i], half[i]), c(-100, i-.5), col="green", lty=2)
        lines(c(-100, half[i]), c(i-.5, i-.5), col="green", lty=2)
        
        # Rasch-Andrich tresholds
        lines(c(rasch[i], rasch[i]), c(-100, rasch_sc[i]), col="blue", lty=2)
        lines(c(-100, rasch[i]), c(rasch_sc[i], rasch_sc[i]), col="blue", lty=2)
      } else {
        # half-tresholds
        lines(c(half[i], half[i]), c(-100, i-.5), col="black", lty=3)
        lines(c(-100, half[i]), c(i-.5, i-.5), col="black", lty=3)
        
        # Rasch-Andrich tresholds
        lines(c(rasch[i], rasch[i]), c(-100, rasch_sc[i]), col="black", lty=4)
        lines(c(-100, rasch[i]), c(rasch_sc[i], rasch_sc[i]), col="black", lty=4)
      }
      
    }
    if(isTRUE(color)) {
      legend("topleft", c("Rasch half-point tresholds", "Rasch-Andrich tresholds"), col=c("green", "blue"), lty=2, inset = .01, bg = "white")
    } else {
      legend("topleft", c("Rasch half-point tresholds", "Rasch-Andrich tresholds"), col="black", lty=c(3,4), inset = .01, bg = "white")
    }
  } else if (show == "info") {
    plot(theta, suminfo, type="l", lwd=3, col="black", ylim=ylim, xlab="trait", ylab="information", main="Information functions")
    for (i in c(1:ncol(pred))) {
      lines(theta, info[,i], lwd=2, col=cols[i])
    }
    legend("topright", legend=c(c(1:ncol(pred))-1, "total"), col=c(cols, "black"), 
           lwd=c(rep(2, ncol(pred)), 3), lty=1, inset=.01, bg = "white")
  } 
  print(rbind(thresholds, half, rasch, rasch_sc))
  invisible(list(pred=pred, thresholds=rbind(thresholds, half, rasch, rasch_sc), info=info, suminfo=suminfo, theta=theta))
}





# Examples pcm.plot ----------------------------------------------------------------

## Example 1
a <- c(1.5, -1.5)
windows(width = 12, height=4.5)
layout(t(c(1:3)))
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = a, show = "trace")
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = a, show = "score")
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = a, show = "info", ylim=c(0,.9))

## example 2, gray scale
b <- c(-1, 1)
windows(width = 12, height=4.5)
layout(t(c(1:3)))
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = b, show = "trace", color=F)
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = b, show = "score", color=F)
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = b, show = "info", color=F)



x <- pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = a, show = "trace")
DescTools::AUC(x$theta, x$suminfo)
y <- pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = b, show = "trace")
DescTools::AUC(y$theta, y$suminfo)




windows(width = 12, height=4.5)
layout(t(c(1:3)))
d <- c(1, -1.3, .3)
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = d, show = "trace")
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = d, show = "score")
pcm.plot(theta=seq(-5,5,by=.001), b=0, thresholds = d, show = "info")



pcm.plot(theta=seq(-5,5,by=.001), b=-1, thresholds = c(-.5, 0, .2, .3), show = "trace", ylim=c(0,.9))
