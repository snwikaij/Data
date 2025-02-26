library(readxl)
library(ggplot2)
library(R2jags)
library(tidyr)
library(MASS)
library(glmmTMB)

#Get wd, this will be the location where all figures will be stored.
wd <- getwd()

#Detect number of cores
numberofcores <- parallel::detectCores()

#Upload the data (literature) and priors
url         <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025_Data_and_priors.xlsx"
destfile    <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
literature  <- read_xlsx(destfile, 1)

#Total dataset with estimates
tot_est <- data.frame(g=literature$Response,
                       s=literature$Type,
                       level=paste(literature$Parameter, literature$Type, literature$Link, literature$Response, sep = "_"),
                       y=literature$estimate,
                       V=literature$estimate_se,
                       m=literature$Link,
                       d=literature$Source,
                       st=literature$DOI)

################################################
#Data sources and systematic search (figure S1)#
################################################

errdat <- data.frame(y=c(0.071,0.500,0.214,0.143,0.071,0.071,0.214,0.071,0.071,0.214,
                         0.571,0.357,0.143,0.143,0.286,0.000,0.286,0.214,0.286,0.143,
                         0.000,0.357,0.286))

me <- summary(glmmTMB(y~1, data=errdat, zi=~1, family = beta_family(link="logit")))
mu <- me$coefficients$cond[1, 1]
mu <- plogis(mu)

va  <- sd(errdat$y)^2
phi <- (mu * (1 - mu) / va) - 1

a <- mu * phi
b <- (1 - mu) * phi

x <- seq(0.0001, 0.9999, 0.0001)
y <- dbeta(x, a, b)

grid.est <- data.frame(x=x,
                       likelihood=y/max(y),
                       prior=(sn::dsn(x, 0.05, 0.5, 12))/max(sn::dsn(x, 0.05, 0.5, 12)))

grid.est$posterior <- grid.est$likelihood*grid.est$prior/max(grid.est$likelihood*grid.est$prior)
post.samp          <- sample(grid.est$x, 10000, T, grid.est$posterior)

fig1               <- tidyr::gather(grid.est,  "key", "value", -x)
cri                <- round(quantile(post.samp, c(.05, .95)),2)

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

cri <- round(hdi_fun(post.samp, .9),2)

p  <- paste0("Prior: mean=", round(1/(1+2),2), " sd=", round(sqrt((1 * 2)/((1 + 2)^2*(1 + 2 + 1))),2))
l  <- paste0("Likelihood: mean=", round(mu,2), " sd=", round(sqrt(va),2))
po <- paste0("Posterior: mean=", round(mean(post.samp),2), " sd=", round(sd(post.samp),2), "\n(5%=", cri[1], ", 95%=", cri[2],")")

fig1$key           <- plyr::mapvalues(fig1$key, from=c("prior", "likelihood", "posterior"),
                                      to=c(p, l, po))

fig1$key           <- factor(fig1$key, levels = c(p,l,po))

figs1 <- ggplot(fig1, aes(x=x, y=value, group=key, col=key))+
  geom_line(lwd=1.2)+ylab("Scaled probability density")+
  xlab("False Negatives")+
  theme_classic()+labs(col="")+
  theme()

ggsave(figs1, filename=paste0(wd,"/Fig_S1.jpeg"), units = "mm", width = 160, height = 80, dpi = 300)

##############################################
#Accumulating estimates estimates (figure S2)#
##############################################

beta        <- tot_est[literature$Parameter == "b1",]
nestimate   <- aggregate(beta, y~g+s+m, length)
nestimate$y <- paste0("n=", nestimate$y)
beta$g      <- factor(beta$g, levels=rev(c("Bacteria", "Algae", "Macrophytes", "Invertebrates", "Fish")))
beta$s      <- factor(beta$s, levels=c("Salinity", "Oxygen", "Sediment", "Thermal", "Flow", "Nutrient-N", "Nutrient-P"))

splitdf     <- split(beta, beta$m)

logdf       <- splitdf$log
logitdf     <- splitdf$logit

plbox1 <- ggplot(logdf, aes(g, y)) +
  coord_flip() +
  geom_hline(yintercept = 0, lty=2, lwd=.2, col="tomato3") +
  geom_boxplot(outlier.shape = NA, alpha=0.8, lwd=0.2) +
  ylim(-2, 2) +
  ylab("Parameter estimate [=regression coefficient]") +
  facet_wrap(.~factor(s,levels=c("Salinity", "Oxygen", "Sediment", "Thermal", "Flow", "Nutrient-N", "Nutrient-P")), scales = "free") +
  geom_point(size=.1, position = position_jitter(width = 0.2), show.legend = F)+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=4),
        strip.text = element_text(size = 6),
        strip.background = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth =0.2),
        axis.ticks = element_line(linewidth =0.2),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom",
        legend.title = element_text(size=6),
        legend.text = element_text(size = 4))

plbox2 <- ggplot(logitdf, aes(g, y)) +
  coord_flip() +
  geom_hline(yintercept = 0, lty=2, lwd=.2, col="tomato3") +
  geom_boxplot(outlier.shape = NA, alpha=0.8, lwd=0.2) +
  ylim(-2, 2) +
  ylab("Parameter estimate [=regression coefficient]") +
  facet_wrap(.~factor(s,levels=c("Salinity", "Oxygen", "Sediment", "Thermal", "Flow", "Nutrient-N", "Nutrient-P")), scales = "free") +
  geom_point(size=.1, position = position_jitter(width = 0.2), show.legend = F)+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=4),
        axis.text = element_text(size=4),
        strip.text = element_text(size = 6),
        strip.background = element_rect(linewidth = 0.4),
        axis.line = element_line(linewidth =0.2),
        axis.ticks = element_line(linewidth =0.2),
        legend.key.size = unit(1, "lines"),
        legend.position = "bottom",
        legend.title = element_text(size=6),
        legend.text = element_text(size = 4))

figs2 <- cowplot::plot_grid(plbox1, plbox2, ncol=1,labels = "auto")

ggsave(figs2, filename=paste0(wd,"/Fig_S2.jpeg"), units = "mm", width = 120, height = 180, dpi = 300)

##########################################################################################################
#A priori bias assessment and quality control (Regression-Based Bias Detection (Egger’s Test)  figure S3)#
##########################################################################################################

x  <- seq(-0.5, 0.5, 0.001)
y1 <- dnorm(x, -0.2, 0.1)
y2 <- dnorm(x, 0.2, 0.1)
y3 <- dnorm(x, 0, 0.025)
y1 <- y1/max(y1)
y2 <- y2/max(y2)
y3 <- y3/max(y3)

dfe <- data.frame(x=c(x,x,x),
                  y=c(y1,y2,y3),
                  t=c(rep("M1a", length(x)),rep("M1b", length(x)), rep("M0", length(x))))

figs3 <- ggplot(dfe, aes(x, y, col=t))+
  geom_line(lwd=1.2)+xlab(expression("Intercept " * beta[0]))+
  scale_color_manual(breaks=c("M1a", "M1b", "M0"),
                     values=c("tomato3", "tomato3", "dodgerblue3"))+
  theme_classic()+ylab("Probability Density")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

ggsave(figs3, filename=paste0(wd,"/Fig_S3.jpeg"), units = "mm", width = 200, height = 80, dpi = 300)

##########################################################################################################
#A priori bias assessment and quality control (Regression-Based Bias Detection (Egger’s Test)  figure S4)#
##########################################################################################################

dfreg <- data.frame(est=literature$estimate[literature$Parameter=="b1"],
                    inv_est=literature$estimate[literature$Parameter=="b1"]/literature$estimate_se[literature$Parameter=="b1"],
                    inv_se=1/literature$estimate_se[literature$Parameter=="b1"],
                    n=literature$n[literature$Parameter=="b1"],
                    inv_n=1/literature$n[literature$Parameter=="b1"],
                    s = factor(literature$Source[literature$Parameter=="b1"],levels=c("Figure", "Table", "Dataset")),
                    r=as.factor(paste(literature$DOI[literature$Parameter=="b1"], literature$`Exact response`[literature$Parameter=="b1"])))

egger <- function(){

  for (i in 1:n){
    y[i]              ~ dnorm(mu[i], t)

    mu[i]            <- intercept[i] + slope[sr[i]] * x[i]
    intercept[i]     <- (b0a[sr[i], r[i]]*m0+b0b[sr[i], r[i]]*m1)}

  for(i in 1:nr){
    slope[i]       ~ dnorm(0, 1/1^2)}

  for(i in 1:ns){
    b00[i]    ~ dnorm(0, 1/.025^2)

    p1[i] <- equals(bd, 1)
    p2[i] <- equals(bd, 2)
    p3[i] <- (-0.2*p1[i]+0.2*p2[i])

    b01[i]    ~ dnorm(p3[i], 1/.1^2)

    for(j in 1:nr){
      b0a[i,j]    ~ dnorm(b00[i], 1/.1^2)
      b0b[i,j]    ~ dnorm(b01[i], 1/.1^2)
    }}

  s     ~ dunif(0, 10)
  t     <- 1/s^2

  bd ~ dcat(c(0.5,0.5))

  d     ~ dcat(pw[1:2])
  m0    <- equals(d, 1)
  m1    <- equals(d, 2)}

reglist2 <- list(x=dfreg$inv_se, y=dfreg$inv_est,  n=nrow(dfreg),
                 r=as.factor(dfreg$r), nr=length(unique(dfreg$r)),
                 sr=as.factor(dfreg$s), ns=length(unique(dfreg$s)),
                 pw=c(0.5, 0.5))

egger_model <-jags.parallel(data = reglist2,
                            model.file = egger,
                            parameters.to.save =  c("b0a", "b0b", "d", "s", "slope"),
                            n.chains = 20,
                            n.thin = 5,
                            jags.seed = 666,
                            n.iter = 3000,
                            n.burnin = 1000)

mcmc      <- egger_model$BUGSoutput$sims.list

ext_b0 <- function(mcmc){

  b0ar <- sapply(1:reglist2$ns, function(x) mcmc$b0a[,x,][mcmc$d == 1,])
  b0br <- sapply(1:reglist2$ns, function(x) mcmc$b0b[,x,][mcmc$d == 2,])
  b0   <- rbind(b0ar, b0br)

  set_names <- levels(reglist2$sr)

  diff_list <- list()

  for (i in 1:(ncol(b0)-1)) {
    for (j in (i+1):ncol(b0)) {
      diff_list[[paste0(set_names[i]," vs. ",set_names[j])]] <- b0[,i]-b0[,j]
    }
  }

  colnames(b0) <- levels(reglist2$sr)
  list(chain=gather(as.data.frame(b0)), comparisons=data.frame(comparison=rep(names(diff_list), each=nrow(b0)), value=as.numeric(unlist(diff_list))))}

b0dens <- ext_b0(mcmc)
b01    <- b0dens$chain

labelegger  <- paste0("\nβ0 = ", round(mean(b01$value),2), ", se = ", round(sd(b01$value),2),
                      "\nBF M1/M0 = ",ifelse(((table(mcmc$d)[1]/table(mcmc$d)[2])/(0.5/0.5))<.01|is.na(table(mcmc$d)[2]), "<.01", round((table(mcmc$d)[1]/table(mcmc$d)[2])/(0.5/0.5),2)),
                      "\nn = ", reglist2$n)

plegger <- ggplot(dfreg, aes(inv_est, inv_se, col=s))+xlim(-10,10)+ylim(0,50)+
  labs(alpha=NA)+
  geom_vline(xintercept = mean(b01$value), lty=1, col="tomato3", lwd=0.6)+
  geom_vline(xintercept = 0, lty=2, col="dodgerblue4", lwd=0.6)+
  ylab("Precision \n1/se")+
  xlab("Parameter estimate/se")+
  labs(alpha = NULL)+
  scale_alpha(guid ="none", range = c(0.9, 0.1))+
  geom_point(data=dfreg, aes(size=log(n), alpha=ifelse(abs(inv_est)>1.96, 0, 1)),  pch=19)+
  scale_size_continuous(guide = "none", range=c(0, 5))+
  annotate("text", x=-7, y=40, label=labelegger, size=3)+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.ticks = element_line(linewidth =0.2),
        legend.key.size = unit(1, "lines"),
        axis.title = element_text(size=6),
        legend.position = "bottom")

b0dens$chain$key  <- factor(b0dens$chain$key, levels=c("Figure", "Table", "Dataset"))
pldensegger <- ggplot(b0dens$chain, aes(value, group=key, fill=key))+
  geom_density(alpha=0.25)+xlab(expression("Estimate " * beta[0]))+ylab("Posterior density")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size=6),
        axis.title.x = element_text(size=6),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

b0dens$comparisons$comparison  <- factor(b0dens$comparisons$comparison,
                                         levels=unique(b0dens$comparisons$comparison))

tabdiff <- aggregate(b0dens$comparisons, value~comparison, function(x) round(c(mean(x), sd(x)),2))
tabdiff <- setNames(cbind.data.frame(tabdiff$comparison, tabdiff$value), c("comparision", "mean", "se"))

pldensegger <- pldensegger+annotation_custom(tableGrob(tabdiff,rows = NULL, theme = ttheme_default(base_size = 6,
                                                                                                   padding = unit(c(1.2, 1.2), "mm"))),
                                             xmax=-0.4, ymax=4)

figs4 <- cowplot::plot_grid(plegger, pldensegger,
                               labels = "auto",
                               ncol=1, rel_heights = c(0.65, 0.35))

ggsave(figs4, filename=paste0(wd,"/Fig_S4.jpeg"), units = "mm", width = 120, height = 140, dpi = 300)

#########################################################################################
#A priori bias assessment and quality control (Z-Value Distribution Analysis  figure S5)#
#########################################################################################

beta        <- tot_est[literature$Parameter == "b1",]
histo       <- data.frame(z=beta$y/beta$V, s=beta$d, v=beta$s, g=beta$g)
histo$v     <- as.factor(histo$v)
zv          <- seq(-10, 10, 0.01)

zvmea       <- aggregate(data=histo, z~v, function(x) mean(x[abs(x)<5]))
zvmin       <- aggregate(data=histo, z~v, function(x) sd(x[abs(x)<5]))
zdf         <- cbind(zvmea, zvmin[-1])

curve_list <- list()
for(i in 1:nrow(zdf)){
  sub_df      <- data.frame(v=zdf[i,1], x=zv, z=dnorm(zv, zdf[i,2], zdf[i,3]))
  sub_df$z    <- sub_df$z/max(sub_df$z)
  curve_list[[i]] <- sub_df}

curve_df <- do.call(rbind, curve_list)
histo$g <- factor(histo$g, levels = c("Bacteria", "Algae", "Macrophytes", "Invertebrates", "Fish"))

plhist <- ggplot(histo, aes(z, fill=g))+xlim(-10, 10)+
  facet_wrap(.~v, scales = "free", ncol=2)+labs(fill="Organism group")+xlab("z-value")+
  ylab("Count")+
  geom_vline(xintercept = c(-1.96, 1.96), col="tomato3", lty=3, lwd=0.8)+
  geom_histogram(col="black",
                 alpha=0.4,
                 binwidth = .5,
                 boundary = 0, closed = "left")+
  theme_classic()

extwd       <- ggplot_build(plhist)$data[[2]]
nnorm       <- merge(curve_df, data.frame(v=levels(histo$v), h=aggregate(data=extwd, y~PANEL, max)[,2]), by="v")
nnorm$z     <- nnorm$z*nnorm$h

figs5 <- plhist+geom_line(data=nnorm, aes(x=x, y=z),lty=2, lwd=1.1,  inherit.aes = F)

ggsave(figs5, filename=paste0(wd,"/Fig_S5.jpeg"), units = "mm", width = 200, height = 240, dpi = 300)

###############################
#Prior formulation (figure S6)#
###############################

#Data from Ruaro et al. (2020)
tab1 <- read_xlsx(destfile, 4)
r1   <- tab1$Fz

#Data from Sabater et al. (2018) Problematic because no sample size is given
#Therefore arbitrary assume 20 as sample size
tab2 <- read_xlsx(destfile, 5)
r2   <- tab2$Fz

#Data from Jackson et al. (2016)
tab3 <- read_xlsx(destfile, 6)
r3   <- tab3$Fz

#Data from Mack et al. (2022)
tab4 <- read_xlsx(destfile, 7)
r4   <- tab4$Fz

#Data from Hering et al. (2006)
tab5 <- read_xlsx(destfile, 8)
r5   <- tab5$Fz

#Data from own exploration
fishz <- read_xlsx(destfile, 9)

#SD for all priors
sdlit <- mean(c(sd(r1),
                sd(r2),
                sd(r3,na.rm = T),
                sd(r4),
                sd(r5),
                sd(fishz$B1)))

##################
#Figure prior set#
##################

x <- seq(-2, 2,0.01)
negative_set <- data.frame(x=x,
                           Prior1=dnorm(x, 0, 10),
                           Prior2=dnorm(x, 0, 0.5),
                           Prior3=dnorm(x, -0.2, 0.1),
                           Prior4=dnorm(x, -0.3, 0.3))

neutral_set <- data.frame(x=x,
                          Prior1=dnorm(x, 0, 10),
                          Prior2=dnorm(x, 0, 0.5),
                          Prior3=dnorm(x, -0.3, 0.3),
                          Prior4=dnorm(x, 0.3, 0.3))

positive_set <- data.frame(x=x,
                           Prior1=dnorm(x, 0, 10),
                           Prior2=dnorm(x, 0, 0.5),
                           Prior3=dnorm(x, 0.2, 0.1),
                           Prior4=dnorm(x, 0.3, 0.3))

negative_set[ , -1] <- apply(negative_set[ , -1], 2, function(col) col / max(col))
neutral_set[ , -1]  <- apply(neutral_set[ , -1], 2, function(col) col / max(col))
positive_set[ , -1] <- apply(positive_set[ , -1], 2, function(col) col / max(col))
negative_set        <- gather(negative_set, "key", "value", -x)
neutral_set         <- gather(neutral_set, "key", "value", -x)
positive_set        <- gather(positive_set, "key", "value", -x)

negative_set$type  <- "Negative set"
neutral_set$type   <- "Neutral set"
positive_set$type  <- "Positive set"

priordir <- rbind(negative_set, neutral_set, positive_set)

plot1 <- ggplot(priordir, aes(x, value, col=key))+
  labs(col=NULL)+facet_wrap(.~type, ncol=3)+
  geom_line(lwd=0.6)+xlab(expression("Regression coefficient " * beta[1]))+
  scale_color_manual(breaks=c("Prior1", "Prior2",  "Prior3", "Prior4"),
                     labels=c("Prior 1", "Prior 2", "Prior 3",  "Prior 4"),
                     values = c( "grey70", "black","orange2", "skyblue3"))+
  ylab("Probability density")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#########################
#Directional explanation#
#########################

order_p <- c("Salinity \nincrease", "Oxygen \ndepletion", "Sediment \nenrichment", "Warming",
             "Flow \ncessation", "N-increase", "P-increase")
order_g <- c("Bacteria", "Algae", "Macrophytes", "Invertebrates", "Fish")

heatmap       <- expand.grid(order_p, order_g)
heatmap$Prior <- as.factor(c(0,0,0,0,0,0,0,#Bacteria
                             -1,0,0,0,0,1,1,#Algae
                             -1,0,0,0,1,-1,-1,#Macrophytes
                             -1,-1,-1,-1,-1,0,0,#Invertebrates
                             -1,-1,-1,1,0,0,0))#Fish

heatmap$Var2 <- factor(heatmap$Var2, levels=rev(order_g))

plot2 <- ggplot(heatmap , aes(x=Var1, y=Var2, fill=Prior)) +
  ggplot2::geom_tile()+
  ggplot2::scale_fill_manual(values = c("tomato3", "grey90", "dodgerblue3"),
                             labels=c("Negative", "Neutral", "Positive")) +
  theme_classic()+labs(fill="Prior set")+
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        axis.text.x.top = element_text(),
        axis.ticks.x.bottom = element_blank())

plot1 <- cowplot::plot_grid(NA, plot1, rel_widths = c(0.05, 0.95))
figs6 <- cowplot::plot_grid(plot1, plot2, ncol=1, labels=c("auto"), rel_heights = c(0.35, 0.65))

ggsave(figs6, filename=paste0(wd,"/Fig_S6.jpeg"), units = "mm", width = 140, height = 100, dpi = 300)
