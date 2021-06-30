#DAPC script
#Jamie Hudson
#15 Feb 2021

library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(extrafont)
library(tidyr)

# Input data --------------------------------------------------------------

crob_genepop <- read.genepop("./data/crob_neutral.gen", ncode = 3)

crob_genepop$pop <- factor(c(rep("BUS", 9), rep("EL", 16), rep("FK", 15),
                               rep("HB", 15), rep("KNY", 15), rep("MEL", 15),
                               rep("NEL", 17), rep("PE", 16), rep("PLY", 15),
                               rep("PO", 16), rep("RAV", 6), rep("SB", 12),
                               rep("TB", 16), rep("TG", 7)))

# xvalDapc to calculate number of PCs to retain ---------------------------


set.seed(999)
system.time(crobx <- xvalDapc(tab(crob_genepop, NA.method= "mean"), pop(crob_genepop))) # 40

set.seed(999)
system.time(crobx <- xvalDapc(tab(crob_genepop, NA.method= "mean"), pop(crob_genepop),
                               n.pca = 30:50, n.rep = 1000,
                               parallel = "multicore", ncpus =4L)) # 33

#View which PCs has lowest MSE and highest mean success
crobx[2:6]


# DAPC with a priori ------------------------------------------------------

dapc.crob <- dapc(crob_genepop, n.pca=33, n.da=5)
pdf(file = "./results/190421_dapc.crob_neutral_colour.pdf",
    width = 4,
    height = 4)
scatter(dapc.crob, scree.da=FALSE, col = c("#cc99ff","#E9C166",
                                           "#cc99ff","#E9BC65","#E9B162","#E8AB60", "#E8A65F", "#E89B5C", "#E8955B", 
                                           "#cc99ff", "#E88A58", "#E77F55", "#E77452", 
                                           "#cc99ff"),
        xax=1, yax=2, cex=2, pch=20, cell=1.5, clab = 1.2)
dev.off()
100*dapc.crob$eig[2]/sum(dapc.crob$eig)

pdf(file = "./results/190421_dapc.crob_neutral_1_3.pdf",
    width = 4,
    height = 4)
scatter(dapc.crob, scree.da=FALSE, col = c("#cc99ff","#E9C166",
                                           "#cc99ff","#E9BC65","#E9B162","#E8AB60", "#E8A65F", "#E89B5C", "#E8955B", 
                                           "#cc99ff", "#E88A58", "#E77F55", "#E77452", 
                                           "#cc99ff"),
        xax=1, yax=3, cex=2, pch=20, cell=1.5, clab = 1.2)
dev.off()


# DAPC without a priori grouping ------------------------------------------

crob_genepop$pop <- factor(c(rep("BUS", 9), rep("EL", 16), rep("FK", 15),
                             rep("HB", 15), rep("KNY", 15), rep("MEL", 15),
                             rep("NEL", 17), rep("PE", 16), rep("PLY", 15),
                             rep("PO", 16), rep("RAV", 6), rep("SB", 12),
                             rep("TB", 16), rep("TG", 7)), c("FK", "BUS", "PO", "TG",
                                                             "NEL", "MEL", "KNY", "PE", 
                                                             "EL", "SB", "TB", "HB", 
                                                             "RAV", "PLY"))

grp <- find.clusters(crob_genepop) #300 pcs and 2 clusters


dapc_crob <- dapc(crob_genepop, grp$grp, n.pca=33, n.da = 1)

scatter(dapc_crob, scree.da=FALSE, col = c("#440D54","#FAE724"), legend = T)

# Grunwald compoplot 

dapc.results <- as.data.frame(dapc_crob$posterior)
dapc.results$pop <- pop(crob_genepop)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

pdf(file = "./results/110521_dapc.crob_compo.pdf",
    width = 10.81,
    height = 2.32,
    family = "Arial")
(compo <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
        geom_bar(stat='identity') + 
        scale_fill_manual(values = c("#FAE724", "#440D54", "blue", "red")) +
        facet_grid(~Original_Pop, 
                   scales = "free", 
                   switch= "x") +
        labs(y = "Posterior membership \nprobability",
             x = "Sample site",
             fill = "Assigned cluster") +
        theme(text = element_text(family = "Arial"),
              axis.text.x = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 10),
              axis.ticks.x = element_blank(),
        )
)
dev.off()
