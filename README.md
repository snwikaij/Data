# Data and supplementary information<br />
<p align="justify">
Welcome to my data repository! I hope this data can help you in satisfying your curiosities or strengthen your article(s). Below in the "read.me" each article is linked with its DOI to the online version of the Journal. The files are sorted in the repository by (1) the name of the Journal, (2) first author and (3) year of publication (Journal_Author_Year). The links/source  are ordered from oldest to newest, per year of publication. Every data source, if I am officially allowed to share, can be accessed via the "Link to Data". Equally, for the Supplementary material via "Link to Supplementary material". The data is in the .xlsx or .csv format, meaning it can be directly loaded into R for simplicity. For example, the csv to "Salinity tolerance of aquatic plants ..." can be loaded into R as: 

```
read.csv(url("https://raw.githubusercontent.com/snwikaij/Data/main/Aquatic_Botany_Kaijser_et_al._2019.csv"), header = T)
```

# 2025
### [Global scale quantification of stressor responses in five riverine organism groups](...)<br />
These data include two scripts:<br />
Main Script – Generates the figures for the main text, including Figures S7, S8, S9, S10, S11 and S12 but excluding Figures 5 and 6. <br />
Appendix script - Produces the figures presented in the appendix. <br />
<br />
Additionally, two Excel documents (Link Data_and_priors and Link to Posterior_results) are included. These contain the necessary data and priors to run both the main and appendix scripts, as well as the acces the posterior estimates from the full analysis.
Most functions used in the scripts are part of the  [EcoPostView](https://snwikaij.github.io/EcoPostView/EcoPostView.html) in R. This R package can be installed via GitHub using the `devtools` package and to use all functions JAGS needs to be installed from https://sourceforge.net/projects/mcmc-jags/ <br /> 
<br />
[Link to Main script](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Main_script.R)<br />
[Link to Appendix script](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Appendix_script.R)<br />
<br />
[Link to Data and priors](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Supplementary_Information_2.xlsx)<br />
[Link to Posterior results](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Supplementary_Information_3.xlsx)<br />
<br />
[Link to Response 1](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Simulates_Precision_Error_Response_1.R)<br />
[Link to Response 2](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_Stochastic_Simulates_Response_2.xlsx)<br />

### [Macrophytes and their sedimentary phosphorus niche in lowland rivers](...) <br />
The data used in this article originates from an exploratory sampling protocol in NRW at the moment it is not published (yet) but it is given in the data availability statement I provide the data online. I assume it will be published somewhere 2025(/2026). I provide five files the full matrix, community matrix (0/1), environmental variable matrix, data used for the GLMM and priors. <br />
<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025.xlsx)<br />
[Link to TP niche code](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._niche_TP.R)<br />
[Link to CCA code](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._2025_cca.R)<br />
[Link to GLMM code](https://github.com/snwikaij/Data/blob/main/Unknown_Kaijser_et_al._GLMM_TP.R)<br />

# 2024
### [Differential associations of five riverine organism groups with multiple stressors](https://www.sciencedirect.com/science/article/pii/S0048969724032522?via%3Dihub) <br />
The data used in this article originates from sampling designed within the RESIST-project. More variables are given in this dataset than used in the publication. Perhaps this might lend itself usefull for other analysisis. As I am only disclosed ot share what I used, therefore only the table used for the PCA and the similairity formated data for the different GLMMs are provided. Not the actuall raw data utilized in further studies. <br />
<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/STOTEN_Kaijser_et_al._2024.xlsx)<br />

# 2022
### [Reach hydromorphology: a crucial environmental variable for the occurrence of riverine macrophytes](https://link.springer.com/article/10.1007/s10750-022-04983-w)<br />
The data used in this article originates from multiple water managing authorities within the state of North Rhine-Westphalia (Germany). Some unclear notations and naming in German is present.<br />
<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/Hydrobiologia_Kaijser_et_al._2022.csv)<br />
[Link to Supplementary information](https://github.com/snwikaij/Data/blob/main/Hydrobiologia_Kaijser_et_al._Supplementary_information.docx)<br />
<br />

### [Environmental ranges discriminating between macrophytes groups in European rivers](https://doi.org/10.1371/journal.pone.0269744)<br />
The data used in this article originates from multiple sources and sometimes it is not exactly clear which laboratory procedure and equipment were used to obtain the data. This is one of the dissadvantages of larger datasets from multiple sources, which need to be accepted and considered. Nonethless, even due to its coarseness, it contains valuable data with relative consistent patterns fitting with others. The sources are mentioned in the supplementary information.<br />
<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/PLOS_One_Kaijser_et_al._2022.csv)<br />
[Link to Supplementary information](https://github.com/snwikaij/Data/blob/main/PLOS_One_Kaijser_et_al._2022_Supplementary_information.docx)<br />
<br />

# 2021
### [The interplay of nutrients, dissolved inorganic carbon and algae in determining macrophyte occurrences in rivers](https://doi.org/10.1016/j.scitotenv.2021.146728)<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/STOTEN_Kaijser_et_al._2021_macrophytes.csv)<br />
[Link to Supplementary information](https://github.com/snwikaij/Data/blob/main/STOTEN_Kaijser_et_al._2021_Supplementary_information.docx)<br />
<br />

# 2019
### [Salinity tolerance of aquatic plants indicated by monitoring data from the Netherlands](https://doi.org/10.1016/j.aquabot.2019.103129)<br />
The data used in this article originates from a waterboard in the Netherlands. I never spend the effort to gather the information of this dataset, the link to the waterboard is https://scheldestromen.nl/. The people responsible for gathering the information are mentioned in the acknowledgement section.<br />
<br />
[Link to Data](https://github.com/snwikaij/Data/blob/main/Aquatic_Botany_Kaijser_et_al._2019.csv)<br />
[Link to Supplementary information](https://github.com/snwikaij/Data/blob/main/Aquatic_Botany_Kaijser_et_al._2019_Supplementary_information.docx)<br />
<br />
</p>






