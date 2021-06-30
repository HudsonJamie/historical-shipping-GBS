library(tidyverse)
library(reshape)
library(corrplot)
library(ggpubr)
library(patchwork)
library(ggtext)


# Microcosmus squamiger ---------------------------------------------------

msqua_table <- read.table("./data/msqua_table_fst.txt")
msqua_table <- msqua_table %>% 
  select(-nw.pacific)
msqua_table <- msqua_table[!rownames(msqua_table) %in% c("nw.pacific"), ]

# 1800

# read in shipping data
ship_network_1800 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1800.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_1800_1 <- ship_network_1800[!rownames(ship_network_1800) %in% c("NW Pacific"), !colnames(ship_network_1800) %in% c("NW.Pacific")]

# 1850 

# read in shipping data
ship_network_1850 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1850.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_1850_1 <- ship_network_1850[!rownames(ship_network_1850) %in% c("NW Pacific"), !colnames(ship_network_1850) %in% c("NW.Pacific")]

# 1900 

# read in shipping data
ship_network_1900 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1900.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_1900_1 <- ship_network_1900[!rownames(ship_network_1900) %in% c("NW Pacific"), !colnames(ship_network_1900) %in% c("NW.Pacific")]

# 1950 

# read in shipping data
ship_network_1950 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1950.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_1950_1 <- ship_network_1950[!rownames(ship_network_1950) %in% c("NW Pacific"), !colnames(ship_network_1950) %in% c("NW.Pacific")]

# 2000 

# read in shipping data
ship_network_2000 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_2000.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_2000_1 <- ship_network_2000[!rownames(ship_network_2000) %in% c("NW Pacific"), !colnames(ship_network_2000) %in% c("NW.Pacific")]

# Tidy the Fst table to match the sites in the shipping data

msqua_matrix <- as.matrix(msqua_table)

msqua_matrix <- round(msqua_matrix,4)
msqua_matrix.melt <- melt(msqua_matrix, na.rm =TRUE)

# Average shipping over time

ship_network_1850_2 <- ship_network_1850_1 %>%  as.data.frame(.) %>% 
  add_column(Aust = 0, .before = 1) %>% 
  add_column(NW.Pacific = 0, .before = 4) %>% 
  add_row(Aust = 0, Med.Sea = 0, NE.Atlantic = 0, NW.Pacific = 0, S.Africa = 0, .before = 1) %>% 
  add_row(Aust = 0, Med.Sea = 0, NE.Atlantic = 0, NW.Pacific = 0, S.Africa = 0, .before = 4)


rownames(ship_network_1850_2)[1] <- "Aust"
rownames(ship_network_1850_2)[4] <- "NE Pacific"

ship_network_1850_2 <- as.matrix(ship_network_1850_2)

x = list(ship_network_1800_1, ship_network_1850_2, ship_network_1900_1, ship_network_1950_1, ship_network_2000_1)
ship_network_ave = reduce(x, `+`) / length(x)


ship_network_ave_1 <- round(ship_network_ave,4)
ship_network_ave.melt <- melt(ship_network_ave_1, na.rm =TRUE)

cbind_all <- cbind(msqua_matrix.melt, ship_network_ave.melt)
colnames(cbind_all) <- c("Var1", "Var2", "Fst", "Var3", "Var4", "ship")

cbind_all_unique <- cbind_all %>% arrange(Fst) %>% 
  group_by(Fst) %>% 
  filter(Var1 != Var2) %>% 
  mutate(ave_ship = mean(ship)) %>% 
  ungroup() %>% 
  mutate(id = rep(1:10, each=2)) %>% 
  distinct(id, .keep_all = TRUE)

shapiro.test(cbind_all_unique$Fst) # W = 0.75544, p-value = 0.000201
shapiro.test(cbind_all_unique$ship) # W = 0.39991, p-value = 4.452e-08

ggqqplot(cbind_all_unique$Fst, ylab = "Fst")
ggqqplot(cbind_all_unique$ship, ylab = "ship")

msqua_all <- ggscatter(cbind_all_unique, x = "ave_ship", y = "Fst", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "spearman",
                       xlab = "Average shipping intensity (1750-2000)", ylab = "Fst",
                       cor.coeff.args = list(label.x = 100, label.sep = "\n"),
                       ggtheme = theme_bw())

# Ciona_robusta -----------------------------------------------------------

crob_table <- read.table("./data/crob_table_fst.txt")

crob_table <- crob_table[!rownames(crob_table) %in% c("ne.pacific"), ]

# 1800 

# read in shipping data
ship_network_1800 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1800.csv", sep=" "))

# remove NE Pacific from shipping data
ship_network_1800_1 <- ship_network_1800[!rownames(ship_network_1800) %in% c("NE Pacific"), !colnames(ship_network_1800) %in% c("NE.Pacific")]

# 1850 

# read in shipping data
ship_network_1850 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1850.csv", sep=" "))

# remove NE Pacific from shipping data
ship_network_1850_1 <- ship_network_1850[!rownames(ship_network_1850) %in% c("NE Pacific"), !colnames(ship_network_1850) %in% c("NE.Pacific")]

# 1900 

# read in shipping data
ship_network_1900 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1900.csv", sep=" "))

# remove NE Pacific from shipping data
ship_network_1900_1 <- ship_network_1900[!rownames(ship_network_1900) %in% c("NE Pacific"), !colnames(ship_network_1900) %in% c("NE.Pacific")]

# 1950 

# read in shipping data
ship_network_1950 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_1950.csv", sep=" "))

# remove NE Pacific from shipping data
ship_network_1950_1 <- ship_network_1950[!rownames(ship_network_1950) %in% c("NE Pacific"), !colnames(ship_network_1950) %in% c("NE.Pacific")]

# 2000 

# read in shipping data
ship_network_2000 <- as.matrix(read.csv("./data/ShipNetwork_TempDev_2000.csv", sep=" "))

# remove NW Pacific from shipping data
ship_network_2000_1 <- ship_network_2000[!rownames(ship_network_2000) %in% c("NE Pacific"), !colnames(ship_network_2000) %in% c("NE.Pacific")]

# Tidy the Fst table to match the sites in the shipping data

crob_matrix <- as.matrix(crob_table)

crob_matrix <- round(crob_matrix,4)
crob_matrix.melt <- melt(crob_matrix, na.rm =TRUE)

# Average shipping over time

ship_network_1850_2 <- ship_network_1850_1 %>%  as.data.frame(.) %>% 
  add_column(Aust = 0, .before = 1) %>% 
  add_row(Aust = 0, Med.Sea = 0, NE.Atlantic = 0, NW.Pacific = 0, S.Africa = 0, .before = 1)  


rownames(ship_network_1850_2)[1] <- "Aust"

ship_network_1850_2 <- as.matrix(ship_network_1850_2)

x = list(ship_network_1800_1, ship_network_1850_2, ship_network_1900_1, ship_network_1950_1, ship_network_2000_1)
ship_network_ave = reduce(x, `+`) / length(x)


ship_network_ave_1 <- round(ship_network_ave,4)
ship_network_ave.melt <- melt(ship_network_ave_1, na.rm =TRUE)

cbind_all <- cbind(crob_matrix.melt, ship_network_ave.melt)
colnames(cbind_all) <- c("Var1", "Var2", "Fst", "Var3", "Var4", "ship")

cbind_all_unique <- cbind_all %>% arrange(Fst) %>% 
  group_by(Fst) %>% 
  filter(Var1 != Var2) %>% 
  mutate(ave_ship = mean(ship)) %>% 
  ungroup() %>% 
  mutate(id = rep(1:10, each=2)) %>% 
  distinct(id, .keep_all = TRUE)

shapiro.test(cbind_all_unique$Fst) # W = 0.75544, p-value = 0.000201
shapiro.test(cbind_all_unique$ship) # W = 0.39991, p-value = 4.452e-08

ggqqplot(cbind_all_unique$Fst, ylab = "Fst")
ggqqplot(cbind_all_unique$ave_ship, ylab = "ship")

crob_all <- ggscatter(cbind_all_unique, x = "ave_ship", y = "Fst", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "spearman",
                      xlab = "Average shipping intensity (1750-2000)", ylab = "Fst",
                      cor.coeff.args = list(label.x = 100, label.sep = "\n"),
                      ggtheme = theme_bw())

msqua_all / crob_all +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.title = ggtext::element_markdown(),
                                axis.text = element_markdown())) 

ggsave(paste0("./figures/both_ave_ship_fst", format(Sys.time(), "%d%m%Y"), ".png"),
       dpi = 320)


