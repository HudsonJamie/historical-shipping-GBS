#hierfstat script to calculate pairwise Fst values of Microcosmus squamiger
#Jamie Hudson
# Created 24 Sep 2019
# Modified 03 March 2021

library(adegenet)
library(hierfstat)
library(reshape2)
library(tidyverse)

# Input data --------------------------------------------------------------

neutral.gen <- read.genepop("../../data/msqua_neutral.gen", ncode = 3) # read in genepop file

neutral.hierfstat <- genind2hierfstat(neutral.gen) #convert genind to hierfstat df

#Rename pop sites so they order correctly in heatmap
neutral.hierfstat$pop <- as.factor(c(rep("A2", 15), rep("A7", 15), rep("A6", 20),
                                     rep("A5", 18), rep("A1", 16), rep("A8", 9),
                                     rep("A9", 13), rep("B1", 16), rep("B2", 15),
                                     rep("C1", 15), rep("B7", 7), rep("A4", 14),
                                     rep("B3", 15), rep("B6", 12), rep("A3", 17),
                                     rep("B9", 10), rep("B4", 13), rep("B8", 13),
                                     rep("C2", 17), rep("B5", 10)))


neutral_pwfst <- genet.dist(neutral.hierfstat, method = "WC84")
str(neutral_pwfst)

neutral_pwfst_matrix<- as.matrix(neutral_pwfst)

colnames(neutral_pwfst_matrix)

write.table(neutral_pwfst_matrix, file="./data/neutral_pwfst_matrix",row.names = T,col.names = T) #rename to something appropriate

# in case you are reading in from file
neutral_pwfst_matrix <- as.matrix(read.table("./data/neutral_pwfst_matrix"))

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
        scale_y_discrete(limits = rev(levels(neutral_pwfst.melt$X2)), labels = rev(c("BU", "AL", "MEL", "MAN", "BF", "AZ", "AR", "CA", "CAD", "CHI", "CU", "MAT", "PB",
                                                                                       "SA", "MB", "KN", "PE", "PA", "EL", "RB"))) +
        scale_x_discrete(labels = c("BU", "AL", "MEL", "MAN", "BF", "AZ", "AR", "CA", "CAD", "CHI", "CU", "MAT", "PB",
                                    "SA", "MB", "KN", "PE", "PA", "EL", "RB")) + 
        ggtitle(expression(italic("Microcosmus squamiger")), paste0(title, "\n", species)) +
        labs( x = "Sampling Site", y = "Sampling Site") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13)) + 
        theme(axis.title = element_text(size = 16),legend.text = element_text(size =20), legend.title = element_text(size =30)) +
        theme(plot.title = element_text(size = 17)) +
        coord_fixed() +
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.title = element_text(size = 14),legend.text = element_text(size =10), legend.title = element_text(size =16)))


# Global Fst --------------------------------------------------------------

msqua_basicstats <- basic.stats(neutral.hierfstat, digits = 4)

msqua_overall <- t(as.data.frame(msqua_basicstats$overall))

# Reg specific Fst --------------------------------------------------------
## Note can't do BF and AZ as they are only one site per region

#rename sites
neutral.hierfstat$pop <- as.factor(c(rep("AL", 15), rep("AR", 15), rep("AZ", 20),
                                 rep("BF", 18), rep("BU", 16), rep("CA", 9),
                                 rep("CAD", 13), rep("CHI", 16), rep("CU", 15),
                                 rep("EL", 15), rep("KN", 7), rep("MAN", 14),
                                 rep("MAT", 15), rep("MB", 12), rep("MEL", 17),
                                 rep("PA", 10), rep("PB", 13), rep("PE", 13),
                                 rep("RB", 17), rep("SA", 10)))


aus <- neutral.hierfstat %>% 
    filter(pop %in% c("AL", "BU", "MEL", "MAN")) %>% 
    mutate(pop = factor(pop))

aus <- basic.stats(aus, digits = 4)
aus_overall <- t(as.data.frame(aus$overall))

introduced <- neutral.hierfstat %>% 
    filter(!pop %in% c("AL", "BU", "MEL", "MAN")) %>% 
    mutate(pop = factor(pop))

introduced <- basic.stats(introduced, digits = 4)
introduced_overall <- t(as.data.frame(introduced$overall))

spain <- neutral.hierfstat %>% 
    filter(pop %in% c("SA", "CAS", "CAD", "CHI", "AR", "MAT", "PB", "CU")) %>% 
    mutate(pop = factor(pop))

spain <- basic.stats(spain, digits = 4)
spain_overall <- t(as.data.frame(spain$overall))

south_africa <- neutral.hierfstat %>% 
    filter(pop %in% c("SB", "TB", "HB", "KNY", "PE", "EL")) %>% 
    mutate(pop = factor(pop))

south_africa <- basic.stats(south_africa, digits = 4)
south_africa_overall <- t(as.data.frame(south_africa$overall))

sessionInfo()
