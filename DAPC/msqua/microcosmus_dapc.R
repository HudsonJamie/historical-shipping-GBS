#DAPC script
#Jamie Hudson
#18 Feb 2021

library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# Input data --------------------------------------------------------------

msqua_genepop <- read.genepop("./data/msqua_neutral.gen", ncode = 3)

msqua_genepop$pop <- factor(c(rep("AL", 15), rep("AR", 15), rep("AZ", 20),
                                rep("BF", 18), rep("BU", 16), rep("CA", 9),
                                rep("CAD", 13), rep("CHI", 16), rep("CU", 15),
                                rep("EL", 15), rep("KNY", 7), rep("MAN", 14),
                                rep("MAT", 15), rep("MB", 12), rep("MEL", 17),
                                rep("PA", 10), rep("PB", 13), rep("PE", 13),
                                rep("RB", 17), rep("SA", 10)))

# xvalDapc to calculate number of PCs to retain ---------------------------


set.seed(999)
system.time(msqua.x <- xvalDapc(tab(msqua_genepop, NA.method= "mean"), pop(msqua_genepop))) # 60

set.seed(999)
system.time(msqua.x <- xvalDapc(tab(msqua_genepop, NA.method= "mean"), pop(msqua_genepop),
                               n.pca = 50:70, n.rep = 1000,
                               parallel = "multicore", ncpus =4L)) # 51

#View which PCs has lowest MSE and highest mean success
msqua.x[2:6]


# DAPC with a priori ------------------------------------------------------

dapc.micro <- dapc(msqua_genepop, n.pca=51, n.da=7)
pdf(file = "microcosmus/results/190421_dapc.msqua_neutral_colours_labels.pdf",
    width = 4,
    height = 4)
scatter(dapc.micro, scree.da=FALSE, col = c("#2a9d8f","#E9C166","#E9BC65","#E9B663",
                                            "#2a9d8f", "#E9B162","#E8AB60","#E8A65F","#E8A05E", "#E89B5C","#E8955B",
                                            "#2a9d8f","#E89059", "#E88A58",
                                            "#2a9d8f","#E88557","#E77F55","#E77A54", "#E77452", "#E76F51"),
        xax=1, yax=2, cex=2, pch=20, cell=1.5, clab = 1.2)
dev.off()

100*dapc.micro$eig[2]/sum(dapc.micro$eig)

pdf(file = "microcosmus/results/190421_dapc.msqua_neutral_1_3.pdf",
    width = 4,
    height = 4)
scatter(dapc.micro, scree.da=FALSE, col = c("#2a9d8f","#E9C166","#E9BC65","#E9B663",
                                            "#2a9d8f", "#E9B162","#E8AB60","#E8A65F","#E8A05E", "#E89B5C","#E8955B",
                                            "#2a9d8f","#E89059", "#E88A58",
                                            "#2a9d8f","#E88557","#E77F55","#E77A54", "#E77452", "#E76F51"),
        xax=1, yax=3, cex=2, pch=20, cell=1.5, clab = 1.2)
dev.off()
100*dapc.micro$eig[2]/sum(dapc.micro$eig)

# DAPC without a priori ---------------------------------------------------

msqua_genepop$pop <- factor(c(rep("AL", 15), rep("AR", 15), rep("AZ", 20),
                              rep("BF", 18), rep("BU", 16), rep("CA", 9),
                              rep("CAD", 13), rep("CHI", 16), rep("CU", 15),
                              rep("EL", 15), rep("KNY", 7), rep("MAN", 14),
                              rep("MAT", 15), rep("MB", 12), rep("MEL", 17),
                              rep("PA", 10), rep("PB", 13), rep("PE", 13),
                              rep("RB", 17), rep("SA", 10)), levels=c("BU", "AL", "MEL", "MAN",
                                                                      "BF", "AZ", "SA", "CA",
                                                                      "CAD", "CHI", "CU",
                                                                      "PB", "MAT", "AR", "MB",
                                                                      "KNY", "PE", "PA",
                                                                      "EL", "RB"))

grp.msqua <- find.clusters(msqua_genepop, n.pca=300) #300 pcs and 2 clusters
dapc_msqua <- dapc(msqua_genepop, grp.msqua$grp, n.pca=51)

pdf(file = "microcosmus/results/dapc.msqua_denovo.pdf",
    width = 4,
    height = 4)
scatter(dapc_msqua, scree.da=FALSE, col = c("orange", "#2a9d8f"), legend = T)
dev.off()

100*dapc_msqua$eig[1]/sum(dapc_msqua$eig)
#Population grouping

table(pop(msqua_genepop))
table(pop(msqua_genepop),grp.msqua$grp)


# Grunwald compoplot ------------------------------------------------------

dapc.results <- as.data.frame(dapc_msqua$posterior)
dapc.results$pop <- pop(msqua_genepop)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

pdf(file = "microcosmus/results/dapc.msqua_compo.pdf",
    width = 10.81,
    height = 2.32,
    family = "Arial")
(compo <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
        geom_bar(stat='identity') + 
        scale_fill_manual(values = c("orange", "#2a9d8f")) +
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


