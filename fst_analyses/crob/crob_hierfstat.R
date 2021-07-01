#hierfstat script to calculate pairwise Fst values of Ciona robusta
#Jamie Hudson
# Created 24 Sep 2019
# Modified 03 March 2021

library(adegenet)
library(hierfstat)
library(reshape2)
library(tidyverse)

# Input data --------------------------------------------------------------

neutral.gen <- read.genepop("./data/crob_neutral.gen", ncode = 3) # read in genepop file

neutral.hierfstat <- genind2hierfstat(neutral.gen) #convert genind to hierfstat df

neutral.hierfstat$pop <- as.factor(c(rep("A2", 9), rep("A9", 16), rep("A1", 15),
                                     rep("B3", 15), rep("A7", 15), rep("A6", 15),
                                     rep("A5", 17), rep("A8", 16), rep("B5", 15),
                                     rep("A3", 16), rep("B4", 6), rep("B1", 12),
                                     rep("B2", 16), rep("A4", 7)))

neutral_pwfst <- genet.dist(neutral.hierfstat, method = "WC84")

str(neutral_pwfst)

neutral_pwfst_matrix<- as.matrix(neutral_pwfst)

write.table(neutral_pwfst_matrix, file="./data/crob_neutral_pwfst_matrix",row.names = T,col.names = T) #rename to something appropriate

# in case you are reading in from file
neutral_pwfst_matrix <- as.matrix(read.table("./data/crob_neutral_pwfst_matrix"))

neutral_pwfst_matrix <- round(neutral_pwfst_matrix,4)
neutral_pwfst.melt <- melt(neutral_pwfst_matrix, na.rm =TRUE)
summary(neutral_pwfst.melt$value)

lt <- lower.tri(neutral_pwfst_matrix)
neutral_pwfst_matrix[!lt] <- 0
neutral_pwfst.melt <- melt(neutral_pwfst_matrix, na.rm =TRUE)
neutral_pwfst.melt$value[neutral_pwfst.melt$value == 0] <- NA
title <- paste("Pairwise FST, WC (1984), Neutral loci")
species <- paste(ncol(neutral.hierfstat) - 1, "loci")

(no.ci_plot <- ggplot(data = neutral_pwfst.melt, aes(X2, X1, fill = value))+ geom_tile(color = "white")+ 
        scale_fill_viridis_c(option = "magma",direction = -1, name = expression(F[ST]), na.value = "white", limits = c(-0.0044, 0.31))  +
        scale_y_discrete(limits = rev(levels(neutral_pwfst.melt$X2)), labels = rev(c("FU", "BU", "PO", "TO",
                                                                                       "NEL", "MEL", 
                                                                                       "KN", "PE", "EL",
                                                                                       "SB", "TB", "HB",
                                                                                       "RAV", "PLY"))) +
        scale_x_discrete(labels = c("FU", "BU", "PO", "TO", "NEL", "MEL", "KN", "PE", "EL", "SB", "TB", "HB","RAV", "PLY")) + 
        ggtitle(expression(italic("Ciona robusta")), paste0(title, "\n", species)) +
        labs( x = "Sampling Site", y = "Sampling Site") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13)) + 
        theme(axis.title = element_text(size = 16),legend.text = element_text(size =20), legend.title = element_text(size =30)) +
        theme(plot.title = element_text(size = 17)) +
        coord_fixed() +
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 14),legend.text = element_text(size =10), legend.title = element_text(size =16)))

# Global Fst --------------------------------------------------------------

crob_basicstats <- basic.stats(neutral.hierfstat, digits = 4)

crob_overall <- t(as.data.frame(crob_basicstats$overall))


# Reg specific Fst --------------------------------------------------------

#rename sites
neutral.hierfstat$pop <- as.factor(c(rep("BUS", 9), rep("EL", 16), rep("FK", 15),
                                     rep("HB", 15), rep("KNY", 15), rep("MEL", 15),
                                     rep("NEL", 17), rep("PE", 16), rep("PLY", 15),
                                     rep("PO", 16), rep("RAV", 6), rep("SB", 12),
                                     rep("TB", 16), rep("TG", 7)))


nwp <- neutral.hierfstat %>% 
    filter(pop %in% c("PO", "TG", "FK", "BUS")) %>% 
    mutate(pop = factor(pop))

nwp <- basic.stats(nwp, digits = 4)
nwp_overall <- t(as.data.frame(nwp$overall))

introduced <- neutral.hierfstat %>% 
    filter(!pop %in% c("PO", "TG", "FK", "BUS")) %>% 
    mutate(pop = factor(pop))

introduced <- basic.stats(introduced, digits = 4)
introduced_overall <- t(as.data.frame(introduced$overall))

southafrica <- neutral.hierfstat %>% 
    filter(pop %in% c("SB", "HB", "TB", "KNY", "PE", "EL")) %>% 
    mutate(pop = factor(pop))

southafrica <- basic.stats(southafrica, digits = 4)
southafrica_overall <- t(as.data.frame(southafrica$overall))

wsa <- neutral.hierfstat %>% 
    filter(pop %in% c("SB", "HB", "TB")) %>% 
    mutate(pop = factor(pop))

wsa <- basic.stats(wsa, digits = 4)
wsa_overall <- t(as.data.frame(wsa$overall))

esa <- neutral.hierfstat %>% 
    filter(pop %in% c("KNY", "PE", "EL")) %>% 
    mutate(pop = factor(pop))

esa <- basic.stats(esa, digits = 4)
esa_overall <- t(as.data.frame(esa$overall))

ausnz <- neutral.hierfstat %>% 
    filter(pop %in% c("MEL", "NEL")) %>% 
    mutate(pop = factor(pop))

ausnz <- basic.stats(ausnz, digits = 4)
ausnz_overall <- t(as.data.frame(ausnz$overall))

eu <- neutral.hierfstat %>% 
    filter(pop %in% c("PLY", "RAV")) %>% 
    mutate(pop = factor(pop))

eu <- basic.stats(eu, digits = 4)
eu_overall <- t(as.data.frame(eu$overall))


# Session info ------------------------------------------------------------

sessionInfo()
