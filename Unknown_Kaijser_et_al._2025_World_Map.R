library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(readxl)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(readxl)

###########
#World map#
###########

#Upload the data (literature) and priors
url         <- "https://raw.githubusercontent.com/snwikaij/Data/main/Unknown_Kaijser_et_al._2025_Supplementary_Information_2.xlsx"
destfile    <- tempfile(fileext = ".xlsx")
download.file(url, destfile, mode = "wb")
literature  <- read_xlsx(destfile, 1)

#Organized to filter out names
store   <- unique(paste0(literature$DOI, "$", literature$Country))
store2  <- do.call(rbind.data.frame, strsplit(store, "\\$"))
colnames(store2) <- c("DOI", "country")
store3           <- store2[!duplicated(store2$DOI),]
countries        <- as.data.frame(table(trimws(unlist(strsplit(store3$country, "\\,")))))
studies_data     <- data.frame(country=countries$Var1, n_studies=countries$Freq)
world            <- ne_countries(scale = "medium", returnclass = "sf")
unique(studies_data$country)[!unique(studies_data$country) %in% unique(world$name)]

#change names
world$name <- plyr::mapvalues(world$name,
                              from=c("Côte d'Ivoire", "United States of America", "United Kingdom", "Trinidad and Tobago", "Bosnia and Herz.", "eSwatini"),
                              to=c("Ivory coast", "USA", "UK", "Trinidad", "Bosnia-Herzegovina", "Eswatini"))
unique(studies_data$country)[!unique(studies_data$country) %in% unique(world$name)]

world_data <- merge(world, studies_data, by.x = "name", by.y = "country", all.x = TRUE)

worldplot <- ggplot(data = world_data) +
  geom_sf(aes(fill = n_studies), color = "white", size = 0.1) +
  theme_classic() +
  labs(fill = "Number \nof studies") +
  scale_fill_gradientn(
    colors = c("lightskyblue1", "dodgerblue4"),
    trans = "log",
    na.value = "grey95",
    breaks = c(1, 5, 10, 20, 40),
    labels = c("1", "5", "10", "20", "40")) +
  theme(legend.position = "right",
        axis.text = element_text(size = 10),
        axis.ticks = element_blank()) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 90),
           expand = FALSE) +
  theme(
    plot.margin = margin(1, 1, 1, 1),
    panel.grid.major = element_line(color = "grey80", linetype = "dotted", linewidth = 0.4),
    panel.grid.minor = element_blank())

#Had no patience anymore
panel2 <- literature
panel2 <- panel2[panel2$Parameter == "b1",]
panel2 <- panel2[c("Type", "Response")]
panel2 < -as.data.frame(table(panel2$Type, panel2$Response))
panel2$Var2 <- factor(panel2$Var2, c("Bacteria", "Algae", "Macrophytes", "Invertebrates", "Fish"))

panel2$Var1 <- plyr::mapvalues(panel2$Var1, from=c("Salinity", "Oxygen", "Sediment", "Thermal", "Flow", "Nutrient-N", "Nutrient-P"),
                               to=c("Salinity-increase", "Oxygen-depletion", "Sediment-enrichment", "Warming", "Flow-cessation", "N-increase", "P-increase"))
panel2$Var1 <- factor(panel2$Var1, rev(c("Salinity-increase", "Oxygen-depletion", "Sediment-enrichment", "Warming", "Flow-cessation", "N-increase", "P-increase")))

pl2 <- ggplot(data=panel2, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", col="black")+coord_flip()+theme_classic()+labs(fill="Group")+
  theme(axis.title = element_blank())+ylab("Instances")+
  theme(legend.text = element_text(size=8),
        axis.text.y = element_text(size=8),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(1, "mm"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2))

fig1plot <- cowplot::plot_grid(cowplot::plot_grid(NULL, worldplot, rel_widths = c(0.05, 0.95)), pl4, ncol=1, labels = "auto")

ggsave(fig1plot, filename="C:/Users/admin/OneDrive/Bureaublad/Paper3/Figures/Fig1_main.pdf", units = "mm", width = 180, height = 120, dpi = 1000)
