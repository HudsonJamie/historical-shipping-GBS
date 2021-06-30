# diveRsity script
# Jamie Hudson
# Created: 12 Sep 2019
# Edited: 08 Mar 2021

library(diveRsity)
library(adegenet)
library(tidyverse)

# Fis ---------------------------------------------------------------------

popgen_stats <- divBasic("../adegenet/microcosmus/data/msqua_neutral.gen", outfile = NULL, gp = 3, 
                         bootstraps = 100)

save(popgen_stats, file = "msqua_popgen_stats.RData") # something suitable
load("divBasic_msqua_neutral_diveRsity.RData")
load("~/jah1g20/OneDrive - University of Southampton/Steves_samples/diveRsity/msqua_popgen_stats.RData")

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
  add_column(pops = c("AL", "AR", "AZ", "BF", "BU", "CA",
               "CAD", "CHI", "CU", "EL", "KNY", "MAN",
               "MAT", "MB", "MEL", "PA", "PB", "PE", 
               "RB", "SA")) %>% 
  rename("Fis" = V1,
         "lower" = V2, 
         "upper" = V3)

ciTable$pop_names <- factor(ciTable$pops, 
                            levels=c("BU", "AL", "MEL", "MAN",
                                     "BF", "AZ", "SA", "CA",
                                     "CAD", "CHI", "CU",
                                     "PB", "MAT", "AR", "MB",
                                     "KNY", "PE", "PA",
                                     "EL", "RB"))

ciTable <- ciTable %>% 
  arrange(pop_names) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 16)))

ciTable$group <- factor(ciTable$group, 
                        levels=c("Native", "Introduced"))

# plot results
# plot results
ggplot(ciTable, aes(x = pop_names, y = Fis)) +
  geom_point(size = 6, aes(colour = group)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = group), width = 0, size = 1.5) +
  scale_color_manual(values = c("#2e9e8f", "#e9c166"),
                     labels = c("Native", "Introduced")) +
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

ggsave(paste0("msqua_Fis_values",
              format(Sys.time(), "%m%Y"),
              ".png"),
       dpi = 320,
       type = "cairo-png")

#Saving 12.9 x 8.35 in image

# Extract Ho and He too
Ho <- as.data.frame(popgen_stats$Ho)
He <- as.data.frame(popgen_stats$He)

mean_Ho <- Ho %>% summarise_all(mean, na.rm = T)
colnames(mean_Ho) <- c("AL", "AR", "AZ", "BF", "BU", "CA",
                       "CAD", "CHI", "CU", "EL", "KN", "MAN",
                       "MAT", "MB", "MEL", "PA", "PB", "PE", 
                       "RB", "SA") # same order as input file
mean_Ho <- t(mean_Ho)
mean_Ho <- as_tibble(mean_Ho)

mean_Ho <- mean_Ho %>% 
  rename("Ho" = V1) %>% 
  add_column(pops = c("BU", "AL", "MEL", "MAN",
                      "BF", "AZ", "SA", "CA",
                      "CAD", "CHI", "CU",
                      "PB", "MAT", "AR", "MB",
                      "KNY", "PE", "PA",
                      "EL", "RB"))

mean_Ho$pops <- factor(mean_Ho$pops, 
                       levels = c("BU", "AL", "MEL", "MAN",
                                  "BF", "AZ", "SA", "CA",
                                  "CAD", "CHI", "CU",
                                  "PB", "MAT", "AR", "MB",
                                  "KNY", "PE", "PA",
                                  "EL", "RB"))

mean_Ho <- mean_Ho %>% 
  arrange(pops) 

msqua_mean_Ho <- mean_Ho %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 16)))

(msqua_ho <- ggplot(msqua_mean_Ho, aes(x = pops, y = Ho)) +
    geom_point(size = 4, aes(colour = group)) +
    scale_color_manual(values = c("#2e9e8f", "#e9c166"),
                       labels = c("Native", "Introduced")) +
    labs(x = "Site names",
         y = "Ho") +
    lims(y = c(0.055, 0.085)) +
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

mean_He <- He %>% summarise_all(mean, na.rm = T)
colnames(mean_He) <- c("AL", "AR", "AZ", "BF", "BU", "CA",
                       "CAD", "CHI", "CU", "EL", "KN", "MAN",
                       "MAT", "MB", "MEL", "PA", "PB", "PE", 
                       "RB", "SA") # same order as input file
mean_He <- t(mean_He)
mean_He <- as_tibble(mean_He)

mean_He <- mean_He %>% 
  rename("He" = V1) %>% 
  add_column(pops = c("BU", "AL", "MEL", "MAN",
                      "BF", "AZ", "SA", "CA",
                      "CAD", "CHI", "CU",
                      "PB", "MAT", "AR", "MB",
                      "KNY", "PE", "PA",
                      "EL", "RB"))

mean_He$pops <- factor(mean_He$pops, 
                       levels = c("BU", "AL", "MEL", "MAN",
                                  "BF", "AZ", "SA", "CA",
                                  "CAD", "CHI", "CU",
                                  "PB", "MAT", "AR", "MB",
                                  "KNY", "PE", "PA",
                                  "EL", "RB"))

mean_He <- mean_He %>% 
  arrange(pops) 

msqua_mean_He <- mean_He %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 16)))

(msqua_he <- ggplot(msqua_mean_He, aes(x = pops, y = He)) +
    geom_point(size = 4, aes(colour = group)) +
    scale_color_manual(values = c("#2e9e8f", "#e9c166"),
                       labels = c("Native", "Introduced")) +
    labs(x = "Site names",
         y = "He") +
    lims(y = c(0.09, 0.14)) +
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
# Count number of private alleles

priv_allele.msqua <- readGenepop(infile = "../adegenet/microcosmus/data/msqua_neutral.gen", gp = 3)
priv_allele.msqua$pop_names <- c("AL", "AR", "AZ", "BF", "BU", "CA",
                                "CAD", "CHI", "CU", "EL", "KN", "MAN",
                                "MAT", "MB", "MEL", "PA", "PB", "PE", 
                                "RB", "SA")

priv_allele.msqua$pop <- as.factor(c(rep("AL", 15), rep("AR", 15), rep("AZ", 20),
                                    rep("BF", 18), rep("BU", 16), rep("CA", 9),
                                    rep("CAD", 13), rep("CHI", 16), rep("CU", 15),
                                    rep("EL", 15), rep("KN", 7), rep("MAN", 14),
                                    rep("MAT", 15), rep("MB", 12), rep("MEL", 17),
                                    rep("PA", 10), rep("PB", 13), rep("PE", 13),
                                    rep("RB", 17), rep("SA", 10)))


allele.freq.msqua <- priv_allele.msqua$allele_freq

length(allele.freq.msqua) # number of alleles

msqua.allele.zeros <- lapply(allele.freq.msqua, function(x) which(x==0)) # Tells us which population and allele combination have a frequency of 0

msqua.private.alleles <- lapply(msqua.allele.zeros, function(x) which(length(x) == length(priv_allele.msqua$pop_names) - 1)) # 12 zeros means only one pop has that allele

length(which(msqua.private.alleles==1)) # Number of private alleles

msqua.private.alleles.list <- allele.freq.msqua[which(msqua.private.alleles==1)]

private.alleles.df <- data.frame(matrix(unlist(msqua.private.alleles.list), nrow=length(msqua.private.alleles.list), byrow=TRUE))

colnames(private.alleles.df) <- as.factor(c(rep("AL", 2), rep("AR", 2), rep("AZ", 2),
                                            rep("BF", 2), rep("BU", 2), rep("CA", 2),
                                            rep("CAD", 2), rep("CHI", 2), rep("CU", 2),
                                            rep("EL", 2), rep("KN", 2), rep("MAN", 2),
                                            rep("MAT", 2), rep("MB", 2), rep("MEL", 2),
                                            rep("PA", 2), rep("PB", 2), rep("PE", 2),
                                            rep("RB", 2), rep("SA", 2)))

private.alleles.df_2 <- private.alleles.df[seq(1, ncol(private.alleles.df),2)]

# For each population count number of private alleles

private_allele <- function(dataset, col_name) {
  col_name <- enquo(col_name)
  dataset %>%
    summarise(count = sum(!!col_name != 0 & !!col_name != 1))
}

private_allele(private.alleles.df_2, AL) #15
private_allele(private.alleles.df_2, AR) #7
private_allele(private.alleles.df_2, AZ) #12
private_allele(private.alleles.df_2, BF) #15
private_allele(private.alleles.df_2, BU) #23
private_allele(private.alleles.df_2, CA) #0
private_allele(private.alleles.df_2, CAD) #1
private_allele(private.alleles.df_2, CHI) #9
private_allele(private.alleles.df_2, CU) #3
private_allele(private.alleles.df_2, EL) #7
private_allele(private.alleles.df_2, KN) #3
private_allele(private.alleles.df_2, MAN) #78
private_allele(private.alleles.df_2, MAT) #8
private_allele(private.alleles.df_2, MB) #6
private_allele(private.alleles.df_2, MEL) #9
private_allele(private.alleles.df_2, PA) #8
private_allele(private.alleles.df_2, PB) #8
private_allele(private.alleles.df_2, PE) #2
private_allele(private.alleles.df_2, RB) #8
private_allele(private.alleles.df_2, SA) #2

sessionInfo()

# allelic_richness --------------------------------------------------------

Ar <- as.data.frame(popgen_stats$Allelic_richness)

mean_Ar <- Ar %>% summarise_all(mean, na.rm = T)
colnames(mean_Ar) <- c("AL", "AR", "AZ", "BF", "BU", "CA",
                       "CAD", "CHI", "CU", "EL", "KN", "MAN",
                       "MAT", "MB", "MEL", "PA", "PB", "PE", 
                       "RB", "SA") # same order as input file

mean_Ar <- t(mean_Ar)
mean_Ar <- as_tibble(mean_Ar)

mean_Ar <- mean_Ar %>% 
  rename("Ar" = V1) %>% 
  add_column(pops = c("BU", "AL", "MEL", "MAN",
                      "BF", "AZ", "SA", "CA",
                      "CAD", "CHI", "CU",
                      "PB", "MAT", "AR", "MB",
                      "KNY", "PE", "PA",
                      "EL", "RB"))

mean_Ar$pops <- factor(mean_Ar$pops, 
                       levels=c("BU", "AL", "MEL", "MAN",
                                "BF", "AZ", "SA", "CA",
                                "CAD", "CHI", "CU",
                                "PB", "MAT", "AR", "MB",
                                "KNY", "PE", "PA",
                                "EL", "RB"))

mean_Ar <- mean_Ar %>% 
  arrange(pops) %>% 
  add_column(group = c(rep("Native", 4), rep("Introduced", 16)))

(msqua_ar <- ggplot(mean_Ar, aes(x = pops, y = Ar)) +
  geom_point(size = 4, aes(colour = group)) +
  scale_color_manual(values = c("#2e9e8f", "#e9c166"),
                     labels = c("Native", "Introduced")) +
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

msqua_ar + crob_ar

ggsave(paste0("both_ar_values",
              format(Sys.time(), "%m%Y"),
              ".png"),
       width = 12,
       height = 6,
       dpi = 320,
       type = "cairo-png")


# AMOVA -------------------------------------------------------------------

popgen_stats <- divBasic("../adegenet/microcosmus/data/msqua_neutral.gen", outfile = NULL, gp = 3, 
                         bootstraps = 100)
