library(readxl)
library(R2jags)
library(ggplot2)
library(see)

#download the data
url      <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025.xlsx"
destfile <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
niche    <- read_xlsx(destfile, 4)
priors   <- read_xlsx(destfile, 5)

#set function
univar_niche <- function(env, taxa, prior_mu=NULL, prior_mu_se=NULL, prior_sd=NULL, prior_sd_se=NULL,
                         n_iter=10000, n_burnin=1000, n_thin=1, n_chain=2, n_seed=666,
                         n_intervals=0.95, n_round=1, post_vjust=1.2, data_vjust=6,
                         xlim=NULL, size=3, fill="grey80", col="black", lwd=0.2, alpha=0.7,
                         ylab="Posterior density", xlab="Environmental gradient", ncol=3){

  if(length(env) != length(taxa)) stop("The lenght of the vectors containing environmental and taxon information are not the same.")

  mod_data   <- list(env=env,
                     taxa=as.factor(taxa),
                     NrT=length(unique(taxa)),
                     N=length(taxa),
                     Prior_mu_a=prior_mu^2/prior_mu_se^2,
                     Prior_mu_b=prior_mu/prior_mu_se^2,
                     Prior_sd_a=prior_sd^2/prior_sd_se^2,
                     Prior_sd_b=prior_sd/prior_sd_se^2)

  dist_mod <- function(){

    # Likelihood
    for (i in 1:N){
      env[i]    ~ dgamma(alpha[i], beta[i])
      alpha[i] <- b0[taxa[i]]^2/sigma0[taxa[i]]^2
      beta[i]  <- b0[taxa[i]]/sigma0[taxa[i]]^2}

    #Priors
    for (j in 1:NrT){
      sigma0[j] ~ dgamma(Prior_sd_a[j], Prior_sd_b[j])
      b0[j]     ~ dgamma(Prior_mu_a[j], Prior_mu_b[j])}

  }

  model <- jags.parallel(data = mod_data,
                         model.file = dist_mod,
                         parameters.to.save =  c("b0", "sigma0"),
                         n.iter = n_iter,
                         n.burnin = n_burnin,
                         n.thin = n_thin,
                         n.chains = n_chain,
                         jags.seed = n_seed)

  #Extract chains
  mcmc_df           <- as.data.frame(model$BUGSoutput$sims.matrix)

  #Remove deviance
  mcmc_df           <- mcmc_df[,colnames(mcmc_df) != "deviance"]

  #Split the chains to extract mu and sd
  split_mcmc <- function(mcmc_df, start, end, names){

    post          <- mcmc_df[,start[1]:end[1]]
    colnames(post)<- names
    mcmc_long     <- gather(post)

    return(mcmc_long)}

  #Raw data
  points     <- data.frame(taxa=mod_data$taxa, value=mod_data$env)

  post_mu       <- split_mcmc(mcmc_df, 1, mod_data$NrT, levels(mod_data$taxa))
  post_sd       <- split_mcmc(mcmc_df, mod_data$NrT+1, ncol(mcmc_df), levels(mod_data$taxa))
  posterior     <- setNames(cbind(post_mu, post_sd[,2]), c("taxa", "mu", "sd"))
  post_sim      <- apply(posterior[c(2:3)],2,as.numeric)
  posterior$sim <- apply(post_sim, 1, function(x) rgamma(1, x[1]^2/x[2]^2, x[1]/x[2]^2))

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

  create_tab <- function(x, int=n_intervals, r=n_round, start="Posterior prediction interval:\n"){
    tab <- aggregate(data=x, value~taxa, function(x) c(round(mean(x)), round(sd(x)), round(hdi_fun(x, int))))
    tab <- setNames(cbind.data.frame(tab[,1], tab$value), c("taxa", "mean", "sd", "ll", "ul"))
    tab$label <- paste0(start, " Mean=",tab$mean, " SD=", tab$sd, " ", round((0.5-int/2)*100,r), "%=", tab$ll, " ", round((0.5+int/2)*100,r), "%=", tab$ul)
    return(tab)
  }

  like_table  <- create_tab(points, start = "Data: ")
  post_table  <- create_tab(data.frame(value=posterior$sim, taxa=posterior$taxa))

  if(is.null(xlim)) xlim <- c(min(points$value), max(points$value))

  fig <- ggplot(posterior, aes(x=taxa, y=sim))+
    facet_wrap(.~taxa, scales = "free_y", ncol = ncol)+ coord_flip()+
    geom_violinhalf(alpha=alpha, fill=fill, col=col, lwd=lwd)+
    geom_point(data=points, aes(x=taxa, y=value), shape="l", inherit.aes = F)+
    geom_text(data=post_table, aes(x=taxa, y=mean(xlim), label=label), size=size, vjust=post_vjust)+
    geom_text(data=like_table, aes(x=taxa, y=mean(xlim), label=label), size=size, vjust=data_vjust)+
    ylim(xlim[1], xlim[2])+
    theme_classic()+xlab(ylab)+ylab(xlab)+
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  return(list(Figure=fig,
              Data=like_table[,-6],
              Posterior=post_table[,-6],
              Model=model))

}

#run the function
resultsTP <- univar_niche(env=as.numeric(niche$TP), taxa=as.factor(niche$name),
                          xlab = expression(paste("Pore water TP"," (","μg"~L^-1,")")),
                          xlim = c(0, 10000), n_iter = 10000, n_chain = 6, n_thin=10, size = 2,
                          post_vjust = 1.2, data_vjust = 6,
                          prior_mu = priors$Mu, prior_mu_se = priors$Se,
                          prior_sd = priors$Mu, prior_sd_se = priors$Mu)

#display the figures
resultsTP$Figure
