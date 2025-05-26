library(readxl)
library(vegan)
library(ggplot2)
library(plyr)

#download the data
url         <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025.xlsx"
destfile    <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
com.mat     <- read_xlsx(destfile, 2)
env.mat     <- read_xlsx(destfile, 3)

#Sample size per taxon
n.taxon      <- data.frame(n=colSums(com.mat>0))
n.taxon$text <- ifelse(n.taxon$n>=5, "bold", "plain")
n.taxon$size <- ifelse(n.taxon$n>=5, 4, 3)

#run cca
my.cca <- vegan::cca(com.mat ~ ., data=env.mat, na.action=na.omit)

#total unconstrained variance
rsqcca <- RsquareAdj(my.cca)
labcca <- paste0("R-squared=", round(rsqcca$r.squared,3))

#summary
sumcca <- summary(my.cca)

#explanation by axis
CCA1         <- paste("Axis 1", paste0(paste0("(", paste0(round(as.data.frame(sumcca$cont$importance)[2,1]*100, 1), "%")), ")"))
CCA2         <- paste("Axis 2", paste0(paste0("(", paste0(round(as.data.frame(sumcca$cont$importance)[2,2]*100, 1), "%")), ")"))

#colours for arrows
vectornames <- as.data.frame(sumcca$biplot)
rownames(vectornames) <- plyr::mapvalues(rownames(vectornames), from=c("`P (S)`", "`P (P)`"), to=c("`TP (S)`", "`TP (P)`"))
colornames  <- plyr::mapvalues(rownames(vectornames),
                               from=c("Velocity", "Depth",
                                      "pH", "HCO3",
                                      "`TP (S)`", "`TP (P)`",
                                      "`NH4 (S)`",  "`NH4 (P)`",
                                      "`NO3 (S)`", "`NO3 (P)`",
                                      "`Spatial distance`", "`Distance to source`"),
                               to=c("blue", "blue",
                                    "green", "green",
                                    "green", "brown",
                                    "green", "brown",
                                    "green", "brown",
                                    "orange", "orange"))

#move species names for beter visibility
speciespoints <- as.data.frame(sumcca$species)
speciespoints$CCA2[4] <- speciespoints$CCA2[4]-0.15
speciespoints$CCA2[8] <- speciespoints$CCA2[8]+0.05
speciespoints$CCA2[10] <- speciespoints$CCA2[10]+0.05
speciespoints$CCA2[13] <- speciespoints$CCA2[13]-0.05

#vectors
vectors  <- as.data.frame(sumcca$biplot)
vecnames <- vectors

vecnames$CCA1[5] <- vecnames$CCA1[5]+0.3
vecnames$CCA1[6] <- vecnames$CCA1[6]+0.3
vecnames$CCA1[12] <- vecnames$CCA1[12]+0.3

vecnames$CCA2[5] <- vecnames$CCA2[5]+0.04
vecnames$CCA2[6] <- vecnames$CCA2[6]+0.04
vecnames$CCA2[12] <- vecnames$CCA2[12]+0.04

#vector names
vnames <- plyr::mapvalues(rownames(vecnames), from=c("`P (S)`", "`P (P)`"), to=c("`TP (S)`", "`TP (P)`"))

#plot results
pl3 <- ggplot(speciespoints, aes(CCA1, CCA2))+xlab(CCA1)+ylab(CCA2)+
  geom_vline(xintercept = 0, lty=2)+geom_hline(yintercept = 0, lty=2)+
  geom_segment(data=vectors, aes(x = 0, y = 0, xend = CCA1*2, yend = CCA2*2),
               col=colornames, lwd=0.6, alpha = 0.7,
               arrow = arrow(length = unit(0.35, "cm")))+
  annotate("text", x=ifelse(vecnames$CCA1 > 0, vecnames$CCA1*2+0.15, vecnames$CCA1*2-0.15),
           y=ifelse(vecnames$CCA2 > 0, vecnames$CCA2*2+0.3, vecnames$CCA2*2-0.3),
           label=vnames, col=colornames, size=3, alpha = 0.7)+
  annotate("text", size=n.taxon$size, x=speciespoints$CCA1, y=speciespoints$CCA2,
           fontface=n.taxon$text, label=rownames(speciespoints), col="grey30")+
  xlim(-1.8,2)+
  theme_classic()+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill=NA, size=.8),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title =  element_text(size=10))
