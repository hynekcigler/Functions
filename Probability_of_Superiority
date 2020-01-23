psup <- function(d=NULL, CL=NULL, plot=F, AUC=F) {
  if (is.null(CL)) {
    CL <- pnorm(d/sqrt(2))
    result <- c(CL=CL)
  } else {
    if (CL > 1 | CL < 0) {
      stop("CL has to be within <0, 1> interval.")
    }
    d <- qnorm(CL)*sqrt(2)
    result <- c(d=d)
  }
  if (plot == T) {
    x <- seq(-5,5, by=.01)
    main <- paste0("Cohen d = ",
                   round(d, 2),
                   ", probability of superiority CL = ",
                   round(CL, 2))
    if (AUC == F) {
      cols <- rainbow(2)
      plot(x, dnorm(x), lwd=2, col=cols[1], type="l", 
           xlab = "score", ylab="density", 
           main = main)
      lines(x, dnorm(x-d), lwd=2, col=cols[2])
    } else {
      plot(pnorm(x), pnorm(x, d), type = "l", lwd=2,
           xlab="percentile in the 1st group", ylab = "percentile in the 2nd group",
           main = main)
    }

  }
  cat(paste0("The effect size in Cohen d = ", 
             round(d, 2),
             " means the probability that a randomly selected subject 
from the second sample achieves higher score than randomly selected subject from the first sample 
(probability of superiority) is equal to ",
             round(CL*100, 2),
             " %. \n\n"))
  return(result)
}
