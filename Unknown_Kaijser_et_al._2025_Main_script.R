library(readxl)
library(ggplot2)
library(parallel)
library(usethis)
library(mime)
library(devtools)
library(gridExtra)
library(grid)

#JAGS needs to be installed from https://sourceforge.net/projects/mcmc-jags/
library(R2jags)

#For further information check https://snwikaij.github.io/EcoPostView/EcoPostView.html
#One needs to install the package from GitHub
devtools::install_github("snwikaij/EcoPostView")
library(EcoPostView)

#Get wd, this will be the location where all figures will be stored.
wd <- getwd()

#Detect number of cores
numberofcores <- parallel::detectCores()

#Upload the data (literature) and priors
url         <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025_Supplementary_Information_2.xlsx"
destfile    <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
literature  <- read_xlsx(destfile, 1)
priors      <- read_xlsx(destfile, 3)

#Total number of articles
nrow(litres2 <- literature[!duplicated(literature$DOI),])

#Systematic and non-systematic
table(literature$sys[!duplicated(literature$DOI)])

#Total number of models
nrow(litres1 <- literature[literature$Parameter == "b0",])

#Log and logit-linear models
table(literature$Link[literature$Parameter == "b0"])

#Total log and logit-linear models
table(literature$Link[literature$Parameter == "b1"])

#Create dataset for the model
mod_data <- data.frame(group=literature$Response,
                       predictor=literature$Fignames,
                       level=paste(literature$Parameter, literature$Type, literature$Link, literature$Response, sep = "_"),
                       estimate=literature$estimate,
                       stderr=literature$estimate_se,
                       linkfun=literature$Link)

#Full model
mod <- meta(estimate = mod_data$estimate,
             stderr = mod_data$stderr,
             parameter = do.call(rbind, strsplit(mod_data$level, "_"))[,1],
             predictor = mod_data$predictor,
             link_function = mod_data$linkfun,
             grouping = mod_data$group,
             prior_mu = as.data.frame(priors[c(2,4,6,8)]),
             prior_mu_se = as.data.frame(priors[c(3,5,7,9)]),
             prior_study_var = 5,
             n_iter = 30000,
             n_thin = 30,
             n_chain= numberofcores)

##############
#Select order#
##############

order_p     <- c("Salinity-increase", "Oxygen-depletion", "Sediment-enrichment", "Warming", "Flow-cessation", "N-increase", "P-increase")
order_g     <- c("Bacteria", "Algae", "Macrophytes", "Invertebrates", "Fish")

####################################################
#Posterior residual bias check (figures S7,8 and 9)#
####################################################

residual_check <- rescheck(mod, order_predictor = order_p, order_group = order_g)

ggsave(residual_check$bias_se, filename=paste0(wd,"/Fig_S7.jpeg"), units = "mm", width = 200, height = 120, dpi = 300)
ggsave(residual_check$bias_se_group, filename=paste0(wd,"/Fig_S8.jpeg"), units = "mm", width = 200, height = 120, dpi = 300)
ggsave(residual_check$bias_se_predictor, filename=paste0(wd,"/Fig_S9.jpeg"), units = "mm", width = 200, height = 120, dpi = 300)

###########
#Inversion#
###########

#I am (as main author) the least happy with this. To bring this to a general audience the
#relation of flow and oxygen was inversed. Hence, the idea is that flow increase and oxygen
#increase are not stressors. Therefore they are inversed by multiplying with -1.
#The issue I have is that there exists no inverse of flow velocity (m/s) and oxygen (mg/l).
#For proper interpretation and prediction it does not work well.
mod_inv <- mod
mod_inv$Estimates$b1$estimate[mod_inv$Estimates$b1$predictor %in% c("Flow-cessation", "Oxygen-depletion")] <- -1*mod_inv$Estimates$b1$estimate[mod_inv$Estimates$b1$predictor %in% c("Flow-cessation", "Oxygen-depletion")]

###################################
#Sensitivity analysis (figure S10)#
###################################

#Sensitivity check
mod_sens <- meta(estimate = mod_data$estimate,
                 stderr = mod_data$stderr,
                 parameter = do.call(rbind, strsplit(mod_data$level, "_"))[,1],
                 predictor = mod_data$predictor,
                 link_function = mod_data$linkfun,
                 grouping = mod_data$group,
                 prior_mu = 0,
                 prior_mu_se = 10,
                 n_iter = 30000,
                 n_thin = 30,
                 n_chain= numberofcores)

mod_sens$Estimates$b1$estimate[mod_sens$Estimates$b1$predictor %in% c("Flow-cessation", "Oxygen-depletion")] <- -1*mod_sens$Estimates$b1$estimate[mod_sens$Estimates$b1$predictor %in% c("Flow-cessation", "Oxygen-depletion")]

sensitivity_check <- senscheck(mod_inv, mod_sens, order_predictor = order_p, order_group = order_g)

ggsave(sensitivity_check$posterior_odds, filename=paste0(wd,"/Fig_S10.jpeg"), units = "mm", width = 200, height = 120, dpi = 300)

############################################################
#Visualizing the results (Posterior density plots figure 2)#
############################################################

postdens       <- pdplot(mod_inv, title_size = 2, point_size = .8, line_width = .7, err_bar_lwd = .5,
                         xlab=c("Pooled parameter estimate [=regression coefficient]"),
                         ylab=c("Probability distribution"),
                         order_predictor = order_p, order_group = order_g)

title_pl_log   <- ggplot()+theme_void()+annotate("text", 0, 0, size = 6, label="Taxonomic richness")
title_pl_logit <- ggplot()+theme_void()+annotate("text", 0, 0, size = 6, label="Evenness")

pdp_combined <- cowplot::plot_grid(title_pl_log,
                                   postdens$posterior_density$log,
                                   title_pl_logit,
                                   postdens$posterior_density$logit,
                                   ncol=1,
                                   rel_heights = c(0.025, 0.225, 0.025, 0.225))

ggsave(pdp_combined, filename=paste0(wd,"/Fig_2.jpeg"), units = "mm", width = 190, height = 240, dpi = 300)

###############################################################
#Visualizing the results (Hypothetical outcome plots figure 3)#
###############################################################

hop_pl1 <- hop(mod, group="Bacteria", predictor = "Warming",
               xlab= expression(Temperature ~ "[" ~ degree * C ~ "]"),
               link_function = "logit", ylab="Bacteria evenness",
               xlim=c(0, 3.4), ylim=c(0.7, 1),
               exp_axis = T, round_x_axis = 0,
               nr_hops = 2500, hop_alpha = 0.05,
               hop_lwd = 0.3, xtextsize = 8,
               ytextsize = 8)

hop_pl2 <- hop(mod, group="Algae", predictor = "N-increase",
             xlab= expression(`Nutrient-N`~ "[" * mg ~ L^-1* "]"),
             link_function = "logit", ylab="Algae evenness",
             xlim=c(0, 3.4), ylim=c(0.1,.9),
             exp_axis = T, round_x_axis = 0,
             nr_hops = 2500, hop_alpha = 0.05,
             hop_lwd = 0.3, xtextsize = 8,
             ytextsize = 8)

hop_pl3 <- hop(mod, group="Invertebrates", predictor = "Sediment-enrichment",
             xlab = c("Fine sediment fraction"), ylab="Invertebrate richness",
             link_function = "log", ylim=c(0,80),
             xlim=c(-4.6, 0), exp_axis = T,
             nr_hops = 2500, hop_alpha = 0.05,
             round_x_axis = 2,
             hop_lwd = 0.3, xtextsize = 8,
             ytextsize = 8)

hop_pl4 <- hop(mod, group="Fish", predictor = "Flow-cessation",
            xlab = expression(`Flow velocity` ~ "[" * m ~ s^-1 * "]"),
            link_function = "log", ylab ="Fish richness",
            xlim=c(-4.6, 0.4), ylim=c(0, 25),
            round_x_axis = 2, exp_axis = T,
            nr_hops = 2500, hop_alpha = 0.05,
            hop_lwd = 0.3, xtextsize = 8,
            ytextsize = 8)

fig3_main  <- cowplot::plot_grid(hop_pl1, hop_pl2, hop_pl3, hop_pl4, labels="auto", ncol=2)
ggsave(fig3_main, filename=paste0(wd,"/Fig_3.jpeg"), units = "mm", width = 120, height = 120, dpi = 300)

#############################################################
#Visualizing the results (Partial dependence plots figure 4)#
#############################################################

pdp1 <- hop(mod,
            group="Algae",
            predictor = c("Salinity-increase", "N-increase"),
            xlab= expression(Salinity ~ "[" * µS ~ cm^-1 * "]"),
            ylab= expression(Nutrient-N ~ "[" * mg ~ L^-1 * "]"),
            gradient_title = "Algae\neveness",
            pdp_resolution = 100,
            legend_position = "top",
            link_function = "logit",
            exp_axis = T,
            round_x_axis = 0,
            round_y_axis = 0,
            xlim=c(3.91, 8.5172),
            ylim=c(-1.5, 3.4),
            xtextsize = 8,
            ytextsize = 8)

pdp1 <- pdp1+theme(legend.title = element_text(size=8),
                   legend.key.size = unit(4, 'mm'),
                   legend.text = element_text(size=7))

pdp2 <- hop(mod,
    group="Invertebrates",
    predictor = c("Sediment-enrichment", "Oxygen-depletion"),
    xlab= "Fine sediment fraction",
    ylab= expression(Oxygen ~ "[" * mg ~ L^-1 * "]"),
    gradient_title = "Invertebrate\nrichness",
    pdp_resolution = 100,
    legend_position = "top",
    link_function = "log",
    exp_axis = T,
    round_x_axis = 2,
    round_y_axis = 0,
    xlim=c(-4.61, -0.92),
    ylim=c(1.61, 2.77),
    xtextsize = 8,
    ytextsize = 8)

pdp2 <- pdp2+theme(legend.title = element_text(size=8),
                   legend.key.size = unit(4, 'mm'),
                   legend.text = element_text(size=7))

fig4_main <- cowplot::plot_grid(pdp1, pdp2, ncol = 2, labels = "auto")
ggsave(fig4_main, filename=paste0(wd,"/Fig_4.jpeg"), units = "mm", width = 120, height = 80, dpi = 300)
