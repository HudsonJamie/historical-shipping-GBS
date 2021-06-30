#hierfstat script to calculate pairwise Fst values
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

# Set NAs to 0?
# neutral.hierfstat[is.na(neutral.hierfstat)] <- 0

#Read in pop names to match what is in ms
neutral.hierfstat$pop <- as.factor(c(rep("BUS", 9), rep("EL", 16), rep("FK", 15),
                                     rep("HB", 15), rep("KNY", 15), rep("MEL", 15),
                                     rep("NEL", 17), rep("PE", 16), rep("PLY", 15),
                                     rep("PO", 16), rep("RAV", 6), rep("SB", 12),
                                     rep("TB", 16), rep("TG", 7)))

neutral.hierfstat$pop <- as.factor(c(rep("A2", 9), rep("B5", 16), rep("A1", 15),
                                     rep("B2", 15), rep("B3", 15), rep("A5", 15),
                                     rep("A6", 17), rep("B4", 16), rep("A8", 15),
                                     rep("A3", 16), rep("A7", 6), rep("A9", 12),
                                     rep("B1", 16), rep("A4", 7)))

ordered_pops <- c("FK", "BUS", "PO", "TG",
                  "NEL", "MEL", 
                  "KNY", "PE", "EL",
                  "SB", "TB", "HB",
                  "RAV", "PLY")

neutral.hierfstat <- neutral.hierfstat %>% arrange(match(pop, ordered_pops))


neutral_pwfst <- genet.dist(neutral.hierfstat, method = "WC84")
str(neutral_pwfst)

neutral_pwfst_matrix<- as.matrix(neutral_pwfst)
neutral_pwfst_matrix <- neutral_pwfst_matrix[ordered_pops,ordered_pops]
write.table(neutral_pwfst_matrix, file="crob_neutral_pwfst_matrix",row.names = T,col.names = T) #rename to something appropriate

# in case you are reading in from file
neutral_pwfst_matrix <- as.matrix(read.table("~/jah1g20/OneDrive - University of Southampton/Steves_samples/hierfstat/crob/crob_neutral_pwfst_matrix"))

neutral_pwfst_matrix <- round(neutral_pwfst_matrix,4)
neutral_pwfst.melt <- melt(neutral_pwfst_matrix, na.rm =TRUE)
summary(neutral_pwfst.melt$value)

# Benjamini-Yekutieli FDR correction - a' = a/sum(1/i) where i = number of comparisons
# count number of comparisons for correction
sum(1:length(levels(neutral.hierfstat$pop))-1) # 190 - note not sum(1:20), as despite 20 pops we don't do a pw comparison of each pop on itself
x <- seq(1, sum(1:length(levels(neutral.hierfstat$pop))-1), 1) # 1 to n of number of comparisons
recip.x <- c(1/x)
sum.recip.x <- sum(recip.x)

# Bootstrap for CIs 
pwfst_boot <- boot.ppfst(dat=neutral.hierfstat, quant = c((0.05/sum.recip.x),(1-(0.05/sum.recip.x))), nboot = 1000)
pwfst_boot_lci <- as.matrix(pwfst_boot[[2]])
pwfst_boot_uci <- as.matrix(pwfst_boot[[3]])

melted_pwfst_lci <- melt(pwfst_boot_lci, na.rm =F)
melted_pwfst_uci <- melt(pwfst_boot_uci, na.rm =F)

melted_pwfst_lci <- as.tibble(melted_pwfst_lci)
melted_pwfst_lci <- rename(melted_pwfst_lci, lci = value)

melted_pwfst_uci <- as.tibble(melted_pwfst_uci)
melted_pwfst_uci <- rename(melted_pwfst_uci, uci = value)

fst_ci <- cbind(melted_pwfst_uci, melted_pwfst_lci)
neutral_pwfst.melt.ci <- cbind(fst_ci, neutral_pwfst.melt)
neutral_pwfst.melt.ci <- neutral_pwfst.melt.ci[-c(1:2,4:5)]
neutral_pwfst.melt.ci <- neutral_pwfst.melt.ci[!(is.na(neutral_pwfst.melt.ci$uci) & neutral_pwfst.melt.ci$Var2 != "C2"), ] #Remove for the plot
neutral_pwfst.melt.ci$value[neutral_pwfst.melt.ci$value == 0] <- NA
nrow(neutral_pwfst.melt.ci[which(neutral_pwfst.melt.ci$uci >0 & neutral_pwfst.melt.ci$lci >0),]) # number significant comparisons

# Produce plot
(pwfst.plot <- ggplot(data = neutral_pwfst.melt.ci, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ 
    scale_fill_viridis_c(option = "magma",direction = -1, name = expression(F[ST]), na.value = "white")  +
    scale_y_discrete(limits = rev(levels(neutral_pwfst.melt$Var2)), labels = rev(c("FU", "BU", "PO", "TO",
                                                                                   "MEL", "NEL", "RAV", "PLY",
                                                                                   "SB", "TB", "HB",
                                                                                   "KN", "PE", "EL"))) +
    scale_x_discrete(labels = c("FU", "BU", "PO", "TO", "MEL", "NEL", "RAV", "PLY", "SB", "TB", "HB","KN", "PE", "EL")) + 
    ggtitle(expression(atop("Pairwise FST, WC (1984), Neutral loci", atop(italic("1,140 loci"), "")))) +
    labs(x = "Sampling Site", y = "Sampling Site") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13)) + 
    theme(axis.title = element_text(size = 16),legend.text = element_text(size =20), legend.title = element_text(size =30)) +
    theme(plot.title = element_text(size = 17)) +
    coord_fixed() +
    geom_text(aes(label = ifelse(uci >0 & lci >0, "*","")), size = 10, colour = "white", vjust = 0.75)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14),legend.text = element_text(size =10), legend.title = element_text(size =16)))


lt <- lower.tri(neutral_pwfst_matrix)
neutral_pwfst_matrix[!lt] <- 0
neutral_pwfst.melt <- melt(neutral_pwfst_matrix, na.rm =TRUE)
neutral_pwfst.melt$value[neutral_pwfst.melt$value == 0] <- NA
title <- paste("Pairwise FST, WC (1984), Neutral loci")
species <- paste(ncol(neutral.hierfstat) - 1, "loci")

(no.ci_plot <- ggplot(data = neutral_pwfst.melt, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ 
        scale_fill_viridis_c(option = "magma",direction = -1, name = expression(F[ST]), na.value = "white", limits = c(-0.0044, 0.31))  +
        scale_y_discrete(limits = rev(levels(neutral_pwfst.melt$Var2)), labels = rev(c("FU", "BU", "PO", "TO",
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
