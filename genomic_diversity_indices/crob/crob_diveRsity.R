# diveRsity script
# Jamie Hudson
# Created: 12 Sep 2019
# Edited: 08 Mar 2021

library(diveRsity)
library(adegenet)
library(tidyverse)
library(ggtext)

# Fis ---------------------------------------------------------------------

popgen_stats <- divBasic("../../data/crob_neutral.gen", outfile = NULL, gp = 3, 
                     bootstraps = 100)

save(popgen_stats, file = "../../data/crob_popgen_stats.RData") # something suitable
# load("../../data/crob_popgen_stats.RData")

###explore results stored in variable popgen_stats
names(popgen_stats)

# popgen_stats contains number of objects but interest is in Fis now.
# Fis contains a dataframe for each population with Fis value and 95% CI
# (standard and bias corrected)

#### First extract the relevant information from our results. Note BC_ refers to bias corrected confidence intervals (http://rstudio-pubs-static.s3.amazonaws.com/12475_c81fc78c66324a40b1509ee7ce3a15f6.html)
ciTable <- lapply(popgen_stats$fis, function(x){
  return(c(x["overall", "fis"], x["overall", "BC_lower_CI"],
           x["overall", "BC_upper_CI"]))
})
### convert this into a dataframe
ciTable <- as.data.frame(do.call("rbind", ciTable))

ciTable <- ciTable %>% 
  add_column(pops = c("BUS", "EL", "FK", "HB", "KNY", "MEL",
                      "NEL", "PE", "PLY", "PO", "RAV", "SB",
                      "TB", "TG")) %>% 
  rename("Fis" = V1,
         "lower" = V2, 
         "upper" = V3)


ciTable$pop_names <- factor(ciTable$pops, 
                            levels=c("FK", "BUS", "PO", "TG",
                                     "NEL", "MEL",
                                     "KNY", "PE", "EL",
                                     "SB", "TB", "HB",
                                     "RAV", "PLY"))

ciTable <- ciTable %>% 
  arrange(pop_names) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 10)))

ciTable$group <- factor(ciTable$group, 
                            levels=c("Native", "Introduced"))

# plot results
ggplot(ciTable, aes(x = pop_names, y = Fis)) +
  geom_point(size = 6, aes(colour = group)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = group), width = 0, size = 1.5) +
  scale_color_manual(values = c("#cc99ff", "#e9c166"),
                     labels = c("Putatively native", "Introduced")) +
  labs(x = "Site names",
       y = "F<sub>IS</sub> value") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_markdown(size = 16),
        axis.title.y = element_markdown(size = 16),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9),
        legend.text = element_text(size = 14))

ggsave(paste0("crob_Fis_values",
              format(Sys.time(), "%m%Y"),
              ".png"),
              dpi = 320,
              type = "cairo-png")
#Saving 12.9 x 8.35 in image


# Extract Ho and He too
Ho <- as.data.frame(popgen_stats$Ho)
He <- as.data.frame(popgen_stats$He)

mean_Ho <- Ho %>% summarise_all(mean, na.rm = T)
colnames(mean_Ho) <- c("BU", "EL", "FU", "HB", "KN", "MEL",
                       "NEL", "PE", "PLY", "PO", "RAV", "SB",
                       "TB", "TO") # same order as input file

mean_Ho <- t(mean_Ho)
mean_Ho <- as_tibble(mean_Ho)

mean_Ho <- mean_Ho %>% 
  rename("Ho" = V1) %>% 
  add_column(pops = c("BUS", "EL", "FK", "HB", "KNY", "MEL",
                      "NEL", "PE", "PLY", "PO", "RAV", "SB",
                      "TB", "TG"))

mean_Ho$pops <- factor(mean_Ho$pops, 
                       levels=c("FK", "BUS", "PO", "TG",
                                "MEL", "NEL", "RAV", "PLY",
                                "SB", "TB", "HB",
                                "KNY", "PE", "EL"))

mean_Ho <- mean_Ho %>% 
  arrange(pops) 

crob_mean_Ho <- mean_Ho %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 10)))

crob_mean_Ho$group <- factor(crob_mean_Ho$group, 
                             levels=c("Native", "Introduced"))

(crob_ho <- ggplot(crob_mean_Ho, aes(x = pops, y = Ho)) +
    geom_point(size = 4, aes(colour = group)) +
    scale_color_manual(values = c("#cc99ff", "#e9c166"),
                       labels = c("Putatively native", "Introduced")) +
    labs(x = "Site names",
         y = "Ho") +
    lims(y = c(0.09, 0.145)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.6),
          panel.grid.major.y = element_line(colour = "black",
                                            size = 0.1),
          panel.grid.major.x = element_blank(),
          axis.title = element_markdown(size = 16),
          axis.title.y = element_markdown(size = 16),
          legend.title = element_blank(),
          legend.position = c(0.85, 0.93),
          legend.text = element_text(size = 15)))




mean_He <- He %>% summarise_all(mean, na.rm = T)
colnames(mean_He) <- c("BU", "EL", "FU", "HB", "KN", "MEL",
                       "NEL", "PE", "PLY", "PO", "RAV", "SB",
                       "TB", "TO") # same order as input file
mean_He <- t(mean_He)
mean_He <- as_tibble(mean_He)



mean_He <- mean_He %>% 
  rename("He" = V1) %>% 
  add_column(pops = c("BUS", "EL", "FK", "HB", "KNY", "MEL",
                      "NEL", "PE", "PLY", "PO", "RAV", "SB",
                      "TB", "TG"))

mean_He$pops <- factor(mean_He$pops, 
                              levels=c("FK", "BUS", "PO", "TG",
                                       "NEL", "MEL", 
                                       "KNY", "PE", "EL",
                                       "SB", "TB", "HB",
                                       "RAV", "PLY"))

mean_He <- mean_He %>% 
  arrange(pops) 

crob_mean_He <- mean_He %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 10)))

crob_mean_He$group <- factor(crob_mean_He$group, 
                             levels=c("Native", "Introduced"))

(crob_he <- ggplot(crob_mean_He, aes(x = pops, y = He)) +
    geom_point(size = 4, aes(colour = group)) +
    scale_color_manual(values = c("#cc99ff", "#e9c166"),
                       labels = c("Putatively native", "Introduced")) +
    labs(x = "Site names",
         y = "He") +
    lims(y = c(0.1, 0.25)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.6),
          panel.grid.major.y = element_line(colour = "black",
                                            size = 0.1),
          panel.grid.major.x = element_blank(),
          axis.title = element_markdown(size = 16),
          axis.title.y = element_markdown(size = 16),
          legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 15)))

indices <- cbind(mean_He, mean_Ho) %>% 
  select(-c(pops))

# Count number of private alleles

priv_allele.crob <- readGenepop(infile = "../../data/crob_neutral.gen", gp = 3)
priv_allele.crob$pop_names <- c("BU", "EL", "FU", "HB", "KN", "MEL",
                                 "NEL", "PE", "PLY", "PO", "RAV", "SB",
                                 "TB", "TO")

priv_allele.crob$pop <- as.factor(c(rep("BU", 9), rep("EL", 16), rep("FU", 15),
                                     rep("HB", 15), rep("KN", 15), rep("MEL", 15),
                                     rep("NEL", 17), rep("PE", 16), rep("PLY", 15),
                                     rep("PO", 16), rep("RAV", 6), rep("SB", 12),
                                     rep("TB", 16), rep("TO", 7)))


allele.freq.crob <- priv_allele.crob$allele_freq

length(allele.freq.crob) # number of alleles

crob.allele.zeros <- lapply(allele.freq.crob, function(x) which(x==0)) # Tells us which population and allele combination have a frequency of 0

crob.private.alleles <- lapply(crob.allele.zeros, function(x) which(length(x) == length(priv_allele.crob$pop_names) - 1)) # 12 zeros means only one pop has that allele

length(which(crob.private.alleles==1)) # Number of private alleles

crob.private.alleles.list <- allele.freq.crob[which(crob.private.alleles==1)]

private.alleles.df <- data.frame(matrix(unlist(crob.private.alleles.list), nrow=length(crob.private.alleles.list), byrow=TRUE))

colnames(private.alleles.df) <- as.factor(c(rep("BU", 2), rep("EL", 2), rep("FU", 2),
                                            rep("HB", 2), rep("KN", 2), rep("MEL", 2),
                                            rep("NEL", 2), rep("PE", 2), rep("PLY", 2),
                                            rep("PO", 2), rep("RAV", 2), rep("SB", 2),
                                            rep("TB", 2), rep("TO", 2)))

private.alleles.df_2 <- private.alleles.df[seq(1, ncol(private.alleles.df),2)]

# For each population count number of private alleles

private_allele <- function(dataset, col_name) {
  col_name <- enquo(col_name)
  dataset %>%
    summarise(count = sum(!!col_name != 0 & !!col_name != 1))
}

private_allele(private.alleles.df_2, BU) #4
private_allele(private.alleles.df_2, EL) #54
private_allele(private.alleles.df_2, FU) #19
private_allele(private.alleles.df_2, HB) #11
private_allele(private.alleles.df_2, KN) #35
private_allele(private.alleles.df_2, MEL) #30
private_allele(private.alleles.df_2, NEL) #62
private_allele(private.alleles.df_2, PE) #58
private_allele(private.alleles.df_2, PLY) #21
private_allele(private.alleles.df_2, PO) #17
private_allele(private.alleles.df_2, RAV) #1
private_allele(private.alleles.df_2, SB) #11
private_allele(private.alleles.df_2, TB) #16
private_allele(private.alleles.df_2, TO) #0

sessionInfo()


# allelic_richness --------------------------------------------------------

crob_Ar <- as.data.frame(popgen_stats$Allelic_richness)

crob_mean_Ar <- crob_Ar %>% summarise_all(mean, na.rm = T)
colnames(crob_mean_Ar) <- c("BU", "EL", "FU", "HB", "KN", "MEL",
                       "NEL", "PE", "PLY", "PO", "RAV", "SB",
                       "TB", "TO") # same order as input file

crob_mean_Ar <- t(crob_mean_Ar)
crob_mean_Ar <- as_tibble(crob_mean_Ar)

crob_mean_Ar <- crob_mean_Ar %>% 
  rename("Ar" = V1) %>% 
  add_column(pops = c("BUS", "EL", "FK", "HB", "KNY", "MEL",
                      "NEL", "PE", "PLY", "PO", "RAV", "SB",
                      "TB", "TG"))

crob_mean_Ar$pops <- factor(crob_mean_Ar$pops, 
                       levels=c("FK", "BUS", "PO", "TG",
                                "MEL", "NEL", "RAV", "PLY",
                                "SB", "TB", "HB",
                                "KNY", "PE", "EL"))

crob_mean_Ar <- crob_mean_Ar %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 10)))

crob_mean_Ar$group <- factor(crob_mean_Ar$group, 
                        levels=c("Native", "Introduced"))

(crob_ar <- ggplot(crob_mean_Ar, aes(x = pops, y = Ar)) +
  geom_point(size = 4, aes(colour = group)) +
    scale_color_manual(values = c("#cc99ff", "#e9c166"),
                       labels = c("Putatively native", "Introduced")) +
  labs(x = "Site names",
       y = "Allelic richness") +
  lims(y = c(0.9, 1.5)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45,
                                   vjust = 0.6),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_markdown(size = 16),
        axis.title.y = element_markdown(size = 16),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.1),
        legend.text = element_text(size = 10)))
