#############################################################
## Functions for Stat111
##
## gen_pdf_df / gen_cdf_df (continuous to discrete for graphing)
## plot_kde (kernel density estimator)
## eval_estimator
## compare_estimators
## fs_pois (an example of recursion in fisher scoring)
## mcmc_mh (metropolis hasting algorithm on uniform)
## perm.test (based on difference in means)
## bootstrap
## bootstrap_t_interval
## jackknife
## 
## edward_kim@college.harvard.edu - Jan. 13, 2020
############################################################

gen_pdf_df <- function(dmodel, support) {
  pdf <- sapply(support,dmodel)
  data.frame(support = support, pdf = pdf)
}

gen_cdf_df <- function(pmodel, support) {
  cdf <- sapply(support,pmodel)
  data.frame(support = support, cdf = cdf)
}

plot_kde <- function(data, kernel = "gaussian", bw=10, lab.x, overlay = F, model, dmodel, support=NULL) {
  require(ggplot2)
  df <- data.frame(data = data)
  if(overlay == T) {
    if(is.null(support)) {
      support <- seq(min(data), max(data), length.out=100)
      }
    df.pdf <- gen_pdf_df(dmodel, support)
    ggplot(data=df, aes(x=data)) + 
      geom_density(kernel=kernel, bw=bw, mapping=aes(color="black")) +
      geom_line(data = df.pdf, mapping = aes(x = support, y = pdf, colour="red"))+
      ylab("probability density") + xlab(lab.x) +
      ggtitle(paste("Kernel Density Estimate Plot \nh = ",bw,", ",kernel," kernel", sep="")) +
      labs(color = "Legend") + 
      scale_color_manual(labels = c("KDE", paste(model,"PDF")), values = c("black", "red"))
  } else {
    ggplot(data=df, aes(x=data)) + 
      geom_density(kernel=kernel, bw=bw, mapping=aes(color="black")) +
      ylab("probability density") + xlab(lab.x) +
      ggtitle(paste("Kernel Density Estimate Plot \nh = ",bw,", ",kernel," kernel", sep="")) +
      labs(color = "Legend") + 
      scale_color_manual(labels = c("KDE"), values = c("black"))
  }
}

### Example KDE with Exponential distribution
# df = read.csv("Data/volcano.csv")
# n = nrow(df)
# l.hat = n/100
# date.char <- paste(df$Year, df$Month, df$Day, sep = "-") 
# date <- as.Date(date.char)
# x <- numeric(n)
# x[1] <- as.numeric(date[1] - as.Date("1899-12-31"))
# for (i in 1:(n-1)) {
#   x[i+1] <- as.numeric(date[i+1]-date[i])
# }
# plot_kde(data = x, lab.x = "interarrival time (days)")
# plot_kde(data = x, lab.x = "interarrival time (days)",
#          overlay = T, model = "Expo", dmodel = function(x) {dexp(x, rate = l.hat / 365)})

plot_ecdf <- function(data, lab.x, lab.title = paste("Empirical CDF Plot"), overlay = F, model, pmodel, legend_label = paste(model,"CDF"), support=NULL) {
  require(ggplot2)
  df <- data.frame(data = data)
  if(overlay == T) {
    if(is.null(support)) {
      support <- seq(min(data), max(data), length.out=100)
    }
    df.cdf <- gen_cdf_df(pmodel, support)
    ggplot(data=df,aes(x=data)) + stat_ecdf(mapping = aes(colour="black"), geom ="step") +
      geom_line(data = df.cdf, mapping = aes(x = support, y = cdf, colour="red"))+
      ylab("cumulative probability density") + xlab(lab.x) +
      ggtitle(lab.title) + labs(color = "Legend") +
      scale_color_manual(labels = c("Empirical CDF", legend_label), values = c("black", "red"))
  } else {
    ggplot(data=df,aes(x=data)) + stat_ecdf(mapping = aes(colour="black"), geom ="step") +
      ylab("cumulative probability density") + xlab(lab.x) +
      ggtitle(lab.title) + labs(color = "Legend") +
      scale_color_manual(labels = c("Empirical CDF"), values = c("black"))
  }
}

### Example KDE with Exponential distribution
# df = read.csv("Data/volcano.csv")
# n = nrow(df)
# l.hat = n/100
# date.char <- paste(df$Year, df$Month, df$Day, sep = "-") 
# date <- as.Date(date.char)
# x <- numeric(n)
# x[1] <- as.numeric(date[1] - as.Date("1899-12-31"))
# for (i in 1:(n-1)) {
#   x[i+1] <- as.numeric(date[i+1]-date[i])
# }
# plot_ecdf(data = x, lab.x = "interarrival time (days)")
# plot_ecdf(data = x, lab.x = "interarrival time (days)",
#           overlay = T, model = "Expo", pmodel = function(x) {pexp(x, rate = l.hat / 365)})
# plot_ecdf(data = x, lab.x = "interarrival time (days)",
#           overlay = T, model = "Expo", pmodel = function(x) {pexp(x, rate = l.hat / 365)},
#           legend_label = paste("Expo CDF\nlambda = 0.00805"),
#           lab.title = "Volcano Eruption Interarrival Times Empirical CDF")

eval_estimator <- function(sigma_hat, sigma_star, n_rep, n_data, rmodel, m_data = NULL) {
  if(is.null(m_data)) {m_data <- replicate(n_rep, rmodel(n_data))}
  estimates <- apply(m_data,2,sigma_hat)
  bias = mean(estimates) - sigma_star
  variance = var(estimates)
  mse = bias^2 + variance
  c(bias,variance,mse)
}

### example evaluation with Normal and sigma^2 MLE
# set.seed(111)
# eval_estimator(
#   sigma_hat = function(data) { 1/length(data)*sum((data-mean(data))^2)},
#   sigma_star = 4, n_rep = 10^4, n_data = 10,
#   rmodel = function(n) {rnorm(n, mean=3,sd=2)}
# )
# set.seed(111)
# eval_estimator(
#   sigma_hat = function(data) { 1/(length(data)-1)*sum((data-mean(data))^2)},
#   sigma_star = 4, n_rep = 10^4, n_data = 10,
#   rmodel = function(n) {rnorm(n, mean=3,sd=2)}
# )
## Example with Normal data matrix prepared and input
# adsfadf <- replicate(10^4, rnorm(10))
# eval_estimator(
#   sigma_hat = function(data) { 1/length(data)*sum((data-mean(data))^2)},
#   sigma_star = 0, m_data = adsfadf,
#   rmodel = function(n) {rnorm(n, mean=3,sd=2)}
# )

compare_estimators <- function(sigma_hat1, sigma_hat2, sigma_star, n_rep, n_data, rmodel, eval_crit = eval_estimator) {
  m_data <- replicate(n_rep, rmodel(n_data))
  est1 <- eval_crit(sigma_hat1, sigma_star, m_data = m_data)
  est2 <- eval_crit(sigma_hat1, sigma_star, m_data = m_data)
  output <- matrix(data = c(est1, est2), nrow = 2)
  rownames(output) <- c("sigma_hat1","sigma_hat2")
  colnames(output) <- c("bias","variance","mse")
  output
}
## Example comparison
# compare_estimators(
#   sigma_hat1 = function(data) {1/length(data)*sum((data-mean(data))^2)},
#   sigma_hat2 = function(data) {1/(length(data)-1)*sum((data-mean(data))^2)},
#   sigma_star = 4, n_rep = 10^4, n_data = 10,
#   rmodel = function(n) {rnorm(n, mean=3,sd=2)}
# )

fs_pois <- function(data, MLE, guess) {
  cat(guess,"\n")
  if (abs(guess-MLE) < 0.01) return(guess)
  else return(fs_pois(data, MLE, guess + (mean(data)*exp(-guess))-1))
}
## Example use
# set.seed(111)
# data <- rpois(25,lambda=2)
# fs_pois(data, log(mean(data)), -1)

## Example recursive function: factorial
# rec_factorial <- function(x) {
#   if (x==0) return(1)
#   else return(x*rec_factorial(x-1))
# }

mcmc_mh <- function(s.dist,x0,chain.length=1000) {
  chain <- c(x0,numeric(chain.length))
  current <- x0
  for(i in 1:chain.length+1) {
    proposal <- runif(1)
    accept <- rbinom(1,1,min(c(s.dist(proposal)/s.dist(current),1)))
    chain[i]<- ifelse(accept==1,proposal, current)
    current <- chain[i]
  }
  return(chain)
}
# Example: Beta distribution
# s <- function(x){x^2*(1-x)}
# t <- mcmc_mh(s,x0=0.5)
# plot_ecdf(data = t[100:1001], lab.x = "x",
#          overlay = T, model = "Beta", pmodel = function(x) {pbeta(x,3,2)},
#          legend_label = paste("Beta CDF (3,2)"),
#          lab.title = "MCMC")

perm.test <- function(group1, group0, nrep) {
  data <- c(group1, group0)
  perms <- replicate(nrep,sample(data,length(data),FALSE))
  diffmean <- function(x){mean(x[1:length(group1)])-mean(x[-(1:length(group1))])}
  ts <- apply(perms,2,diffmean)
  return(ts)
}
#### Example application ####
# x <- rnorm(50,0,1)
# y <- rnorm(90,0,10)
# dist <- perm.test(x,y,10^4)
# hist(dist, breaks = 50)
# abline(v = mean(x)-mean(y), col = "blue", lwd = 2)

bootstrap <- function(data, test.statistic, B) {
  boot_samples <- t(replicate(B, sample(data, length(data), replace = TRUE)))
  boot_thetas <- apply(boot_samples, 1, test.statistic) 
  return(c(test.statistic(data),sd(boot_thetas)))
}

bootstrap_t_interval <- function(data, test.statistic, B, alpha=0.05) {
  theta.hat <- test.statistic(data)
  boot.samples <- t(replicate(B, sample(data, length(data), replace = TRUE)))
  theta.tilde <- apply(boot_samples, 1, function(x){bootstrap(x,test.statistic,B)})
  t <- (theta.tilde[1,]-theta.hat)/theta.tilde[2,]
  list.out <- list("CI" = theta.hat - quantile(t,c(1-alpha/2,alpha/2),names=F) * sd(data), "t.dist"= t)
  return(list.out)
}

#### Example ####
# set.seed(111)
# y <- rbeta(30,2,5)
# bootstrap_t_interval(y,median,10^4)

jackknife <- function(data, test.statistic) {
  n <- length(data)
  t <- test.statistic(data)
  jreps <- matrix(numeric(0), nrow=n, ncol=n-1)
  for(i in 1:n) {jreps[i,] <- data[-i]}
  tbar <- 1/n * sum(apply(jreps,1,test.statistic))
  t - (n-1) * (tbar-t)
}
