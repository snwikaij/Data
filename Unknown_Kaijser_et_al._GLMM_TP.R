library(readxl)
library(R2jags)
library(ggplot2)
library(mcmcplots)

#download the data
url      <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025.xlsx"
destfile <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
env      <- read_xlsx(destfile, 3)
full     <- read_xlsx(destfile, 1)

#model
glmm_gamma_jags <- function(){

  # Likelihood
  for (i in 1:N){
    y[i]    ~ dgamma(shape, shape/mu[i])

    log(mu[i])  <- b0 + b1 * x[random[i]]}

  # Random effects
  for (r in 1:n_rand){x[r] ~ dnorm(7.2, 1/7.2^2)}

  # Residuals
  resid<-y-mu

  # Priors
  b0         ~ dnorm(-1, 1/.5^2)
  b1         ~ dnorm(0.69, 1/.14^2)
  shape      ~ dunif(0, 2)}

#model data
data_mod         <- list(N=nrow(env),
                         n_rand=length(unique(full$Location)),
                         random=as.factor(full$Location),
                         x=log(env$`P (P)`),
                         y=env$`P (S)`)

#run model
Jags_mod <- jags.parallel(data = data_mod,
                          model.file = glmm_gamma_jags,
                          parameters.to.save = c("b0", "b1","shape","resid"),
                          jags.seed = 666,
                          n.iter = 20000,
                          n.chains = 8,
                          n.burnin = 1000,
                          n.thin = 20)

#check the chains
mcmcplots::traplot(Jags_mod, parms = c("b1", "b0", "shape"))

#check residuals
residuals <- colMeans(Jags_mod$BUGSoutput$sims.list$resid)
plot(density(residuals))

#pooled estimation for b0
print(b0 <- mean(Jags_mod$BUGSoutput$sims.list$b0))

#pooled estimation for b1
print(b1 <- mean(Jags_mod$BUGSoutput$sims.list$b1))

#create function to extract hdi intervals
hdi_fun              <- function(x, level){

  orddata   <- sort(x)
  nord      <- length(x)
  infomass  <- ceiling(level*nord)
  outmass   <- nord-infomass

  min_width <- Inf
  ll        <- NA
  ul        <- NA

  for(i in 1:(outmass+1)){
    int_width <- orddata[i+infomass-1]-orddata[i]

    if(int_width < min_width){
      min_width <- int_width
      ll <- orddata[i]
      ul <- orddata[i+infomass-1]}}

  c(ll, ul)}

#extract intervals
intb0 <- hdi_fun(Jags_mod$BUGSoutput$sims.list$b0, .95)
intb1 <- hdi_fun(Jags_mod$BUGSoutput$sims.list$b1, .95)

#create label
labgraph  <- paste0(paste0("β0=", round(mean(Jags_mod$BUGSoutput$sims.list$b0),2),
                           " (SE=", round(sd(Jags_mod$BUGSoutput$sims.list$b0),2),
                           "; 2.5%=", round(intb0[1],2),
                           "; 97.5%=", round(intb0[2],2),")"),
                    paste0("\nβ1=", round(mean(Jags_mod$BUGSoutput$sims.list$b1),2),
                           " (SE=", round(sd(Jags_mod$BUGSoutput$sims.list$b1),2),
                           "; 2.5%=", round(intb1[1],2),
                           "; 97.5%=", round(intb1[2],2),")"))

#formulate a df with posterior draws
posterior_draws <-    data.frame(b0=rowMeans(Jags_mod$BUGSoutput$sims.list$b0),
                                 b1=Jags_mod$BUGSoutput$sims.list$b1)

#setup a gradient for hops lines
x_gradient      <- seq(2,9,.05)
hops            <- list()
set.seed(123)
ll              <- sample(1:nrow(posterior_draws), 500)

#generate hops
for(i in ll){

  b0 <- posterior_draws[i,1]
  b1 <- posterior_draws[i,2]

  hops[[i]] <- cbind(i, x_gradient, exp(b0+b1*x_gradient))

}

#place hops in df
hops_realized      <- setNames(do.call(rbind.data.frame, hops), c("l", "x", "y"))
mu                 <- data.frame(x=x_gradient, y=exp(mean(posterior_draws[,1])+mean(posterior_draws[,2])*x_gradient))

#generate figure
ggplot(hops_realized, aes(x, y, group=l))+
  annotate("text", x=4, y=400, size=3, label=labgraph, fontface=2)+
  geom_point(data=data.frame(x=data_mod$x, y=data_mod$y), aes(x, y), shape=19, size=3.8, alpha=0.15, fill="grey20", inherit.aes = F)+
  geom_line(lwd=.6, alpha=0.1, col="dodgerblue2")+xlim(2,9)+ylim(0,550)+
  geom_line(data=mu, aes(x=x, y=y), inherit.aes = F, lwd=1.2, col="dodgerblue4")+
  xlab(expression(paste("Pore water TP [Log transformed]"," (","μg"~L^-1,")")))+
  ylab(expression(paste("Surface water TP "," (","μg"~L^-1,")")))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
