# quick look at m/z vs rt for each lipid class to watch progression with C# and DB#
# looking for outliers
<<<<<<< HEAD

library(ggplot2)
library(tidyverse)



#setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/Nano_Take_Two/First_Third/")

data <- read.csv("Nano_First_Neg_Take_Two_LOBSTAHS_screened_peakdata_2019-06-28T10-46-27_AM-0400.csv")
# get a random column and assign a new column with just
# C number and DB number to label everything
lipidclass <- data %>%
  mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+")))


# going through the big nine plus TAGs and DAGs
fngGSL <- lipidclass %>%
  filter(lipid_class == "fungalGSL")

ggplot(fngGSL, aes(x = peakgroup_rt, y = peakgroup_mz))+
  geom_point()+
  ggtitle("M/Z vs. RT in fngGSL")+
  geom_text(aes(label = compound_name, hjust = 1, vjust = 2))+
  theme(legend.title=element_blank())+
  ylim(650, 900)+
  xlim(500, 1200)

dLCB_GSL_No_FA_OH<- lipidclass %>%
  filter(species == "dLCB_GSL_No_FA_OH", degree_oxidation == 2)

ggplot(dLCB_GSL_No_FA_OH, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in dLCB_GSL_No_FA_OH")+
  theme(legend.title=element_blank())+
  ylim(650, 900)+
  xlim(500, 1200)

TAG <- lipidclass %>%
  filter(species == "TAG")

ggplot(TAG, aes(x = peakgroup_rt, y = peakgroup_mz))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in TAG")+
  theme(legend.title=element_blank())

BLL <- lipidclass %>%
  filter(species == "BLL")

ggplot(BLL, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in BLL")+
  theme(legend.title=element_blank())


DAGunder40 <- lipidclass %>%
 filter(species == "DAG", degree_oxidation == 0, FA_total_no_C < 40)

ggplot(DAGunder40, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag, size = QE004350.mzXML))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in DAGs")+
  theme(legend.title=element_blank())

PC40up <- lipidclass %>%
  filter(species == "PC", degree_oxidation == 0, FA_total_no_C >= 40)

ggplot(PC40up, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in PC40up")+
  theme(legend.title=element_blank())

DGTS_DGTA <- lipidclass %>%
  filter(species == "DGTS_DGTA")

ggplot(DGTS_DGTA, aes(x = peakgroup_rt, y = peakgroup_mz, color = FA_total_no_DB))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in DGTS_DGTA")+
  theme(legend.title=element_blank())

DGCC <- lipidclass %>%
  filter(species == "DGCC")

ggplot(DGCC, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in DGCC")+
  theme(legend.title=element_blank())

LPG <- lipidclass %>%
  filter(species == "LPG")

ggplot(LPG, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag, size = QE004350.mzXML))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in LPG")+
  theme(legend.title=element_blank())

PC <- lipidclass %>%
  filter(species == "PC", degree_oxidation == 0)

ggplot(PC, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in PC")+
  theme(legend.title=element_blank())+
  scale_color_manual(values = c("10%_rtv"="#66CD00", "Red"="#FF3030","ms2v"="#0000FF", "5%_rtv"="#2aff00","Double_Peak?"="#ff9e44","Unknown"="#000000"))

PE <- lipidclass %>%
  filter(species == "PE")

ggplot(PE, aes(x = peakgroup_rt, y = peakgroup_mz, color = FA_total_no_DB))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in PE")+
  theme(legend.title=element_blank())


PG <- lipidclass %>%
  filter(species == "PG")

ggplot(PG, aes(x = peakgroup_rt, y = peakgroup_mz, color = FA_total_no_C))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in PG")+
  theme(legend.title=element_blank())

SQDG <- lipidclass %>%
  filter(species == "SQDG")

ggplot(SQDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in SQDG")+
  theme(legend.title=element_blank())

DGDG <- lipidclass %>%
  filter(species == "DGDG")

ggplot(DGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in DGDG")+
  theme(legend.title=element_blank())

MGDG <- lipidclass %>%
  filter(species == "MGDG")

ggplot(MGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = FA_total_no_DB))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in MGDG")+
  theme(legend.title=element_blank())

MGDG_over_35 <- lipidclass %>%
  filter(species == "MGDG", FA_total_no_C > 35, degree_oxidation < 3)

ggplot(MGDG_over_35, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("KimT M/Z vs. RT in MGDG > 35 C, degree_oxidation < 3")+
  theme(legend.title=element_blank())

MGDG_under_35 <- lipidclass %>%
  filter(species == "MGDG", FA_total_no_C < 35, FA_total_no_C>25, degree_oxidation ==0)

ggplot(MGDG_under_35, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("KimT M/Z vs. RT in MGDG < 35 C")+
  theme(legend.title=element_blank())


TAG_under_50 <- lipidclass %>%
  filter(species == "TAG", FA_total_no_C > 40, FA_total_no_C <= 50)

ggplot(TAG_under_50, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("Nicole M/Z vs. RT in TAG > 40 <= 50C")+
  theme(legend.title=element_blank())

TAG_over_50 <- lipidclass %>%
  filter(species == "TAG", FA_total_no_C > 50, FA_total_no_C < 60)

ggplot(TAG_over_50, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("Nicole M/Z vs. RT in Tag > 50C")+
  theme(legend.title=element_blank())

TAG_under_40 <- lipidclass %>%
  filter(species == "TAG", FA_total_no_C <= 40)

ggplot(TAG_under_40, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("Nicole M/Z vs. RT in Tag <= 40C")+
  theme(legend.title=element_blank())

TAG_over_60 <- lipidclass %>%
  filter(species == "TAG", FA_total_no_C > 60)

ggplot(TAG_over_60, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("Nicole M/Z vs. RT in Tag > 60C")+
  theme(legend.title=element_blank())

# +
#   scale_colour_gradient2(low = "white", mid = "red",
#                          high = "blue", midpoint = 400000000, space = "Lab",
#                          na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
#   xlim(500, 1500)

# used this to compile FFAs, then put them into the csv for quantification
UniqueFFA <- filtered %>% filter(species == "FFA")
OxFFA <- oxidized %>% filter(species == "FFA")
FFAs <- rbind(UniqueFFA, OxFFA)



FFAs <- data %>%
  filter(species == "FFA", degree_oxidation == 0) %>%
  mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+")))


FFAsUnder20 <-  FFAs %>% filter(FA_total_no_C <= 20)
FFAsOver20 <-  FFAs %>% filter(FA_total_no_C > 20)

ggplot(FFAs, aes(peakgroup_rt, peakgroup_mz))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("FFAs m/z vs. rt")+
  theme(legend.title=element_blank())

ggplot(FFAsUnder20, aes(peakgroup_rt, peakgroup_mz, color = QE006421.mzXML))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("FFAs m/z vs. rt")+
  theme(legend.title=element_blank())

ggplot(FFAsOver20, aes(peakgroup_rt, peakgroup_mz, color = QE006421.mzXML))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("FFAs m/z vs. rt")+
  theme(legend.title=element_blank())



UniqueMGDG <- filtered %>% filter(species == "MGDG")
OxMGDG <- oxidized %>% filter(species == "MGDG")
MGDGs <- rbind(UniqueMGDG, OxMGDG)

UniqueMGDG <- UniqueMGDG %>%
  mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))

MGDGs <- MGDGs %>%
  mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))

MGDG30Down <- MGDGs %>% filter(FA_total_no_C < 31)
MGDG31to35 <- MGDGs %>% filter(FA_total_no_C > 31, FA_total_no_C < 36)
MGDG36to40 <- MGDGs %>% filter(FA_total_no_C > 35, FA_total_no_C < 41)
MGDG41to45 <- MGDGs %>% filter(FA_total_no_C > 40, FA_total_no_C < 46)
MGDGover45 <- MGDGs %>% filter(FA_total_no_C > 45)

ggplot(UniqueMGDG, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)

ggplot(MGDG30Down, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)

ggplot(MGDG31to35, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)

ggplot(MGDG36to40, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)

ggplot(MGDG41to45, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)

ggplot(MGDGover45, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
  geom_point()+
  geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
  ggtitle("MGDG m/z vs. rt")+
  theme(legend.title=element_blank())+
  xlim(250, 1500)


## checking out nano marine snow duplicates
MGDG <- dupl_mass %>%
  mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))) %>% filter(species == "MGDG")
ggplot(MGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color =  nLCQE_01340_Jon_Marine_snow_expts._snow_str._15_D8_T6))+
  geom_point()+
  geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
  ggtitle("M/Z vs. RT in MGDG Dupls")+
  theme(legend.title=element_blank())


#############################################
# going through the big lipids plus TAGs, DAGs, and FFAs
# and plotting each

#plotting the rt_dbase
# dbase <- RT_Factor_Dbase %>%
#   filter(is.na(Mean_DNPPE_Factor) == FALSE, is.na(FA_total_no_C) == FALSE)%>%
#   mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+")))




lipidclass <- data %>%
  mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))) %>%
  filter(is.na(FA_total_no_C) == FALSE, degree_oxidation == 0)


lipid_classes <- unique(lipidclass$species)


for (i in 1:length(lipid_classes)){
  Lipid <- lipidclass %>%
    filter(species == paste(lipid_classes[i]))

  print(ggplot(Lipid, aes(x = peakgroup_rt, y = LOBdbase_mz, color =  FA_total_no_C))+
          geom_point()+
          geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
          ggtitle(paste0("M/Z vs. RT in ", lipid_classes[i])))
        # +
        #   scale_color_manual(values = c("10%_rtv"="#66CD00", "Red"="#FF3030","ms2v"="#2aff00", "5%_rtv"="#2aff00","Double_Peak?"="#ff9e44","Unknown"="#000000")))

  ggsave(filename = paste0(lipid_classes[i], "_Nano_First_Third_MZRT.tiff"),
         plot = last_plot(),
         device = "tiff",
         width = 22, height = 17)



}
=======
#
# library(ggplot2)
# library(tidyverse)
#
#
#
# #setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/GSL Tests/NAAMES/")
#
# #data <- read.csv("NAAMES_First_Half_Coded.csv")
# # get a random column and assign a new column with just
# # C number and DB number to label everything
# lipidclass <- data %>%
#   mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))
#
#
# # going through the big nine plus TAGs and DAGs
# fngGSL <- lipidclass %>%
#   filter(lipid_class == "fungalGSL")
#
# ggplot(fngGSL, aes(x = peakgroup_rt, y = peakgroup_mz))+
#   geom_point()+
#   ggtitle("M/Z vs. RT in fngGSL")+
#   geom_text(aes(label = compound_name, hjust = 1, vjust = 2))+
#   theme(legend.title=element_blank())+
#   ylim(650, 900)+
#   xlim(500, 1200)
#
# dLCB_GSL_No_FA_OH<- lipidclass %>%
#   filter(species == "dLCB_GSL_No_FA_OH", degree_oxidation == 2)
#
# ggplot(dLCB_GSL_No_FA_OH, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in dLCB_GSL_No_FA_OH")+
#   theme(legend.title=element_blank())+
#   ylim(650, 900)+
#   xlim(500, 1200)
#
# TAG <- lipidclass %>%
#   filter(species == "TAG")
#
# ggplot(TAG, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in TAG")+
#   theme(legend.title=element_blank())
#
# BLL <- lipidclass %>%
#   filter(species == "BLL")
#
# ggplot(BLL, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in BLL")+
#   theme(legend.title=element_blank())
#
#
# DAGunder40 <- lipidclass %>%
#  filter(species == "DAG", degree_oxidation == 0, FA_total_no_C < 40)
#
# ggplot(DAGunder40, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag, size = QE004350.mzXML))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in DAGs")+
#   theme(legend.title=element_blank())
#
# PC40up <- lipidclass %>%
#   filter(species == "PC", degree_oxidation == 0, FA_total_no_C >= 40)
#
# ggplot(PC40up, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in PC40up")+
#   theme(legend.title=element_blank())
#
# DGTS_DGTA <- lipidclass %>%
#   filter(species == "DGTS_DGTA")
#
# ggplot(DGTS_DGTA, aes(x = peakgroup_rt, y = peakgroup_mz, color = FA_total_no_DB))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in DGTS_DGTA")+
#   theme(legend.title=element_blank())
#
# DGCC <- lipidclass %>%
#   filter(species == "DGCC")
#
# ggplot(DGCC, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in DGCC")+
#   theme(legend.title=element_blank())
#
# LPG <- lipidclass %>%
#   filter(species == "LPG")
#
# ggplot(LPG, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag, size = QE004350.mzXML))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in LPG")+
#   theme(legend.title=element_blank())
#
# PC <- lipidclass %>%
#   filter(species == "PC", degree_oxidation == 0)
#
# ggplot(PC, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in PC")+
#   theme(legend.title=element_blank())+
#   scale_color_manual(values = c("10%_rtv"="#66CD00", "Red"="#FF3030","ms2v"="#0000FF", "5%_rtv"="#2aff00","Double_Peak?"="#ff9e44","Unknown"="#000000"))
#
# PE <- lipidclass %>%
#   filter(species == "PE")
#
# ggplot(PE, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in PE")+
#   theme(legend.title=element_blank())+
#   scale_color_manual(values = c("10%_rtv"="#66CD00", "Red"="#FF3030","ms2v"="#0000FF", "5%_rtv"="#2aff00","Double_Peak?"="#ff9e44","Unknown"="#000000"))
#
#
#
# PG <- lipidclass %>%
#   filter(species == "PG")
#
# ggplot(PG, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in PG")+
#   theme(legend.title=element_blank())
#
# SQDG <- lipidclass %>%
#   filter(species == "SQDG")
#
# ggplot(SQDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in SQDG")+
#   theme(legend.title=element_blank())
#
# DGDG <- lipidclass %>%
#   filter(species == "DGDG")
#
# ggplot(DGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in DGDG")+
#   theme(legend.title=element_blank())
#
# MGDG <- lipidclass %>%
#   filter(species == "MGDG")
#
# ggplot(MGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color = code))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in MGDG")+
#   theme(legend.title=element_blank())
#
# MGDG_over_35 <- lipidclass %>%
#   filter(species == "MGDG", FA_total_no_C > 35, degree_oxidation < 3)
#
# ggplot(MGDG_over_35, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("KimT M/Z vs. RT in MGDG > 35 C, degree_oxidation < 3")+
#   theme(legend.title=element_blank())
#
# MGDG_under_35 <- lipidclass %>%
#   filter(species == "MGDG", FA_total_no_C < 35, FA_total_no_C>25, degree_oxidation < 3)
#
# ggplot(MGDG_under_35, aes(x = peakgroup_rt, y = peakgroup_mz, color = Flag))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("KimT M/Z vs. RT in MGDG < 35 C, degree_oxidation < 3")+
#   theme(legend.title=element_blank())
#
#
# TAG_under_50 <- lipidclass %>%
#   filter(species == "TAG", FA_total_no_C > 40, FA_total_no_C <= 50)
#
# ggplot(TAG_under_50, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("Nicole M/Z vs. RT in TAG > 40 <= 50C")+
#   theme(legend.title=element_blank())
#
# TAG_over_50 <- lipidclass %>%
#   filter(species == "TAG", FA_total_no_C > 50, FA_total_no_C < 60)
#
# ggplot(TAG_over_50, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("Nicole M/Z vs. RT in Tag > 50C")+
#   theme(legend.title=element_blank())
#
# TAG_under_40 <- lipidclass %>%
#   filter(species == "TAG", FA_total_no_C <= 40)
#
# ggplot(TAG_under_40, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("Nicole M/Z vs. RT in Tag <= 40C")+
#   theme(legend.title=element_blank())
#
# TAG_over_60 <- lipidclass %>%
#   filter(species == "TAG", FA_total_no_C > 60)
#
# ggplot(TAG_over_60, aes(x = peakgroup_rt, y = peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("Nicole M/Z vs. RT in Tag > 60C")+
#   theme(legend.title=element_blank())
#
# # +
# #   scale_colour_gradient2(low = "white", mid = "red",
# #                          high = "blue", midpoint = 400000000, space = "Lab",
# #                          na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
# #   xlim(500, 1500)
#
# # used this to compile FFAs, then put them into the csv for quantification
# UniqueFFA <- filtered %>% filter(species == "FFA")
# OxFFA <- oxidized %>% filter(species == "FFA")
# FFAs <- rbind(UniqueFFA, OxFFA)
#
#
#
# FFAs <- coded %>%
#   filter(species == "FFA") %>%
#   mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))
#
#
# FFAsUnder20 <-  FFAs %>% filter(FA_total_no_C <= 20)
# FFAsOver20 <-  FFAs %>% filter(FA_total_no_C > 20)
#
# ggplot(FFAs, aes(peakgroup_rt, peakgroup_mz, color = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("FFAs m/z vs. rt")+
#   theme(legend.title=element_blank())
#
# ggplot(FFAsUnder20, aes(peakgroup_rt, peakgroup_mz, color = QE006421.mzXML))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("FFAs m/z vs. rt")+
#   theme(legend.title=element_blank())
#
# ggplot(FFAsOver20, aes(peakgroup_rt, peakgroup_mz, color = QE006421.mzXML))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("FFAs m/z vs. rt")+
#   theme(legend.title=element_blank())
#
#
#
# UniqueMGDG <- filtered %>% filter(species == "MGDG")
# OxMGDG <- oxidized %>% filter(species == "MGDG")
# MGDGs <- rbind(UniqueMGDG, OxMGDG)
#
# UniqueMGDG <- UniqueMGDG %>%
#   mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))
#
# MGDGs <- MGDGs %>%
#   mutate(C_DB_OX = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"),":", str_extract(degree_oxidation, "\\d+")))
#
# MGDG30Down <- MGDGs %>% filter(FA_total_no_C < 31)
# MGDG31to35 <- MGDGs %>% filter(FA_total_no_C > 31, FA_total_no_C < 36)
# MGDG36to40 <- MGDGs %>% filter(FA_total_no_C > 35, FA_total_no_C < 41)
# MGDG41to45 <- MGDGs %>% filter(FA_total_no_C > 40, FA_total_no_C < 46)
# MGDGover45 <- MGDGs %>% filter(FA_total_no_C > 45)
#
# ggplot(UniqueMGDG, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
# ggplot(MGDG30Down, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
# ggplot(MGDG31to35, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
# ggplot(MGDG36to40, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
# ggplot(MGDG41to45, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
# ggplot(MGDGover45, aes(x = peakgroup_rt, y = peakgroup_mz, fill = degree_oxidation))+
#   geom_point()+
#   geom_text(aes(label = C_DB_OX, hjust = 1, vjust = 2))+
#   ggtitle("MGDG m/z vs. rt")+
#   theme(legend.title=element_blank())+
#   xlim(250, 1500)
#
#
# ## checking out nano marine snow duplicates
# MGDG <- dupl_mass %>%
#   mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))) %>% filter(species == "MGDG")
# ggplot(MGDG, aes(x = peakgroup_rt, y = peakgroup_mz, color =  nLCQE_01340_Jon_Marine_snow_expts._snow_str._15_D8_T6))+
#   geom_point()+
#   geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
#   ggtitle("M/Z vs. RT in MGDG Dupls")+
#   theme(legend.title=element_blank())
#
#
# #############################################
# # going through the big lipids plus TAGs, DAGs, and FFAs
# # and plotting each
# dbase <- RT_Factor_Dbase %>%
#   filter(is.na(Mean_DNPPE_Factor) == FALSE, is.na(FA_total_no_C) == FALSE)%>%
#   mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+")))
#
# lipidclass <- coded %>%
#   mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))) %>%
#   filter(is.na(FA_total_no_C) == FALSE, degree_oxidation == 0)
#
#
# lipid_classes <- unique(lipidclass$species)
#
#
# for (i in 1:length(lipid_classes)){
#   Lipid <- lipidclass %>%
#     filter(species == paste(lipid_classes[i]))
#
#   print(ggplot(Lipid, aes(x = peakgroup_rt, y = LOBdbase_mz, color =  Flag))+
#           geom_point()+
#           geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
#           ggtitle(paste0("M/Z vs. RT in ", lipid_classes[i]))+
#           scale_color_manual(values = c("10%_rtv"="#66CD00", "Red"="#FF3030","ms2v"="#2aff00", "5%_rtv"="#2aff00","Double_Peak?"="#ff9e44","Unknown"="#000000")))
#
#   ggsave(filename = paste0(lipid_classes[i], "_NAAMES_FirstHalf_MZRT.tiff"),
#          plot = last_plot(),
#          device = "tiff",
#          width = 22, height = 17)
#
#
#
# }
>>>>>>> 66722c6bdd1160dfa2676c79ae7f1eae16cdc98a


