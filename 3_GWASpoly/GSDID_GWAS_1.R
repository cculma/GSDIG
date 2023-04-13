library(GWASpoly)
library(ggplot2)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
# set_pheno ---------------------------------------------------------------

pheno <- read.csv("pheno_1500.csv")
trait1 <- colnames(pheno)[2:(length(colnames(pheno))-6)]
summary(pheno$stem)
hist(pheno$color)
summary(pheno$color)
str(pheno)

pheno <- read.csv("Stem_strength_1.csv", row.names = 1)
pheno <- read.csv("Stem_strength.csv", row.names = 1)

# pheno <- read.csv("Stem_strength_1.csv", row.names = 1) # original presentation April 10

head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-5)]
trait1
trait2 <- trait1[-2]

# set_params --------------------------------------------------------------
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")
params <- set.params(fixed = c("PC1","PC2","PC3", "PC4", "PC5"), fixed.type = c("numeric","numeric","numeric","numeric","numeric"), n.PC = 5)

# N <- 490 #Population size
# params <- set.params(geno.freq = 1 - 5/N, fixed = "experiment", fixed.type = "factor")
#
params <- set.params(fixed = c("rep","PC1","PC2","PC3", "PC4", "PC5"), fixed.type = c("factor","numeric","numeric","numeric","numeric","numeric"), n.PC = 5, P3D = F)


data_b <- read.GWASpoly(ploidy=4, pheno.file="pheno_1500.csv",
                        geno.file="DAI21_GWAS_1.txt", format="numeric", n.traits=length(trait1), delim=",")

data_c <- set.K(data_b, LOCO = T)
data_d <- GWASpoly(data_c, models = models_1, traits = "stem", params = params, n.core = 8)

?GWASpoly
?set.K
?read.GWASpoly
?set.params

data_e <- set.threshold(data_d, method= "Bonferroni", level=0.05)
data_e1 <-  get.QTL(data_e)
manhattan.plot(data = data_e)

manhattan.plot(data = data_e, traits = "stage", chrom = 4)

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
write.csv(data_e1, "GWAS1.csv", quote = F, row.names = F)


data_f <- set.threshold(data_d, method= "FDR", level=0.05)
data_f1 <-  get.QTL(data_f)
data_f2 <-  get.QTL(data_f, traits = "stage")
dplyr::count(data_f1, Marker)


# myplot2 <- manhattan.plot(data = data_e) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.y = element_text(size = 12), plot.tag = element_blank(), strip.text.x = element_blank(), strip.background = element_rect(fill = "white", color = "white"), panel.grid.minor = element_blank()) + scale_y_continuous(breaks = seq(0:6))

myplot2 <- manhattan.plot(data = data_e) + theme(legend.position = "none", axis.title.y = element_text(size = 12), plot.tag = element_blank(), strip.background = element_rect(fill = "white", color = "white"), panel.grid.minor = element_blank())

ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/man_1.pdf", plot = myplot2, dpi = 300, width = 6, height = 6, device = cairo_pdf)

dplyr::count(data_e1, Marker, Trait)

# plot LD -----------------------------------------------------------------

p <- LD.plot(data_d)

myplot1 <- p + theme_classic(base_family = "Arial", base_size = 12) + scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 2)) + geom_vline(xintercept = 10, linetype=2, color = "grey") + geom_hline(yintercept = 0.012, linetype=2, color = "grey")

myplot1

ggsave(filename = "~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/4_plots/LD_decay.pdf", plot = myplot1, dpi = 300, width = 4, height = 4, device = cairo_pdf)

# myplot3 <- ggarrange(myplot1, myplot2, nrow = 2, labels = c("a", "b"))


# r^2 ---------------------------------------------------------------------
# phenotypic variance explained (PVE)

data_d
cc <- count(data_e1,Trait)
lev4 <- cc$Trait
lev4
cc1 <- count(data_e1, Model)
cc1$Model

QTL_3 <- data_e1 %>% dplyr::filter(!Model %in% c("diplo-general", "diplo-additive"))


fit_05 <- fit.QTL(data=data_d, trait = "stage",
                  qtl=data_e1[,c("Marker","Model")])

setwd("~/medin297@umn.edu - Google Drive/My Drive/GSDIG-selected/3_GWASpoly/")
write.csv(fit_05, "r2.csv", quote = F, row.names = F)



