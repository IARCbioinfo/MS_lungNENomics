

### MOFA Exp.Meth.Alt.CNV lungNET ###

library(MOFA2)
library(reshape2)
library(ggplot2)
library(cowplot)
library(grid)
library(pracma)
library(plyr)
library(stringr)
library(GGally)

## Load datasets
load("D_expr_red_lungNENomicsCombined.RData")
load("D_meth_red_lungNENomicsCombined.RData")
load("D_Alt_lungNENomicsCombined.RData")
load("D_cnv_seg_lungNENomicsCombined.RData")

load("combined_public_lungNENomics_technical_data.RData")
technical <- technical_complete[which(technical_complete$cohort=="lungNENomicsCombined"),]

## Match sample names
samples <- sort(unique(c(colnames(D_expr_red_273),colnames(D_meth_red_247),colnames(D_Alt),colnames(D_cnv_seg)))) # n=319 
all(samples %in% technical$sample_id) # TRUE
all(technical$sample_id %in% samples) # TRUE
rm(technical_complete)

all(technical$sample_id[which(technical$WGS=="YES")] %in% colnames(D_Alt)) # TRUE
all(colnames(D_Alt) %in% technical$sample_id[which(technical$WGS=="YES")]) # TRUE 
all(technical$sample_id[which(technical$WGS=="YES")] %in% colnames(D_cnv_seg)) # TRUE
all(colnames(D_cnv_seg) %in% technical$sample_id[which(technical$WGS=="YES")]) # TRUE
all(technical$sample_id[which(technical$RNAseq=="YES")] %in% colnames(D_expr_red_273)) # TRUE
all(colnames(D_expr_red_273) %in% technical$sample_id[which(technical$RNAseq=="YES")]) # TRUE
all(technical$sample_id[which(technical$DNAmeth=="YES")] %in% colnames(D_meth_red_247)) # TRUE
all(colnames(D_meth_red_247) %in% technical$sample_id[which(technical$DNAmeth=="YES")]) # TRUE

## Create input datasets with all samples 
D_expr_MOFA <- matrix(NA, nrow(D_expr_red_273), length(which(!samples %in% colnames(D_expr_red_273))), dimnames = list(rownames(D_expr_red_273), samples[which(!samples %in% colnames(D_expr_red_273))]))
D_expr_MOFA <- as.data.frame(D_expr_MOFA) # 5000 x 46
D_expr_MOFA <- merge(D_expr_red_273, D_expr_MOFA, by.x="row.names", by.y="row.names") # 5000x320
D_expr_MOFA <- data.frame(D_expr_MOFA[,-1],row.names=D_expr_MOFA[,1]) # 5000x319
D_expr_MOFA <- D_expr_MOFA[,order(match(colnames(D_expr_MOFA),samples))]
dim(D_expr_MOFA) # 5000 x 319

D_met_MOFA <- matrix(NA, nrow(D_meth_red_247), length(which(!samples %in% colnames(D_meth_red_247))), dimnames = list(rownames(D_meth_red_247), samples[which(!samples %in% colnames(D_meth_red_247))]))
D_met_MOFA <- as.data.frame(D_met_MOFA) # 5000 x 72
D_met_MOFA <- merge(D_meth_red_247, D_met_MOFA, by.x="row.names", by.y="row.names")
D_met_MOFA <- data.frame(D_met_MOFA[,-1], row.names=D_met_MOFA[,1])
D_met_MOFA <- D_met_MOFA[, order(match(colnames(D_met_MOFA), samples))]
dim(D_met_MOFA) # 5000 x 319

D_alt_MOFA <- matrix(NA, nrow(D_Alt), length(which(!samples %in% colnames(D_Alt))), dimnames = list(rownames(D_Alt), samples[which(!samples %in% colnames(D_Alt))]))
D_alt_MOFA <- as.data.frame(D_alt_MOFA) # 163 x 217
D_alt_MOFA <- merge(D_Alt, D_alt_MOFA, by.x="row.names", by.y="row.names")
D_alt_MOFA <- data.frame(D_alt_MOFA[,-1], row.names=D_alt_MOFA[,1])
D_alt_MOFA <- D_alt_MOFA[, order(match(colnames(D_alt_MOFA), samples))]
dim(D_alt_MOFA) # 163 x 319

D_cnv_MOFA <- matrix(NA, nrow(D_cnv_seg), length(which(!samples %in% colnames(D_cnv_seg))), dimnames = list(rownames(D_cnv_seg), samples[which(!samples %in% colnames(D_cnv_seg))]))
D_cnv_MOFA <- as.data.frame(D_cnv_MOFA) # 8 x 217
D_cnv_MOFA <- merge(D_cnv_seg, D_cnv_MOFA, by.x="row.names", by.y="row.names")
D_cnv_MOFA <- data.frame(D_cnv_MOFA[,-1], row.names=D_cnv_MOFA[,1])
D_cnv_MOFA <- D_cnv_MOFA[, order(match(colnames(D_cnv_MOFA), samples))]
dim(D_cnv_MOFA) # 8 x 319

all(colnames(D_expr_MOFA)==samples) # TRUE
all(colnames(D_met_MOFA)==samples) # TRUE
all(colnames(D_alt_MOFA)==samples) # TRUE
all(colnames(D_cnv_MOFA)==samples) # TRUE

save(D_expr_MOFA, D_met_MOFA, D_alt_MOFA, D_cnv_MOFA, file="Inputs_MOFA.RData")

## Create and prepare MOFA object
# Plot density of each input (use the original inputs without NA values)
pdf("Views_Density.pdf",10,7)
par(mfrow=c(2,2),mar = rep(2, 4))
plot(density(as.matrix(D_expr_red_273), na.rm=T, from=min(D_expr_red_273, na.rm=T)), xlab="Gene expression (nrc)", ylab="Frequency", main=paste0("Density on RNA-seq data (n=", ncol(D_expr_red_273), " ; N=", nrow(D_expr_red_273), ")"))
plot(density(as.matrix(D_meth_red_247)), xlab="CpG methylation (M-values)", ylab="Frequency", main=paste0("Density on Meth data (n=", ncol(D_meth_red_247), " ; N=", nrow(D_meth_red_247), ")"))
plot(density(as.matrix(D_Alt), na.rm=T, from=min(D_Alt, na.rm=T)), xlab="Alteration status", ylab="Frequency", main=paste0("Density on SNV/SV data (n=", ncol(D_Alt), " ; N=", nrow(D_Alt), ")"))
plot(density(as.matrix(D_cnv_seg), na.rm=T, from=min(D_cnv_seg, na.rm=T)), xlab="Copy number value", ylab="Frequency", main=paste0("Density on copy number segment data (n=", ncol(D_cnv_seg), " ; N=", nrow(D_cnv_seg), ")"))
dev.off()

D_expr_MOFA  <- as.matrix(D_expr_MOFA)
D_met_MOFA <- as.matrix(D_met_MOFA)
D_alt_MOFA <- as.matrix(D_alt_MOFA)
D_cnv_MOFA <- as.matrix(D_cnv_MOFA)

MOFAobjectB <- create_mofa(list("RNA" = D_expr_MOFA, "Meth" = D_met_MOFA, "Alt" = D_alt_MOFA, "CNV" = D_cnv_MOFA))
data_opts <- get_default_data_options(MOFAobjectB)
model_opts <- get_default_model_options(MOFAobjectB)
model_opts$likelihoods["Alt"] <- "bernoulli"
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobjectB)
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 10000

# Prepare
MOFAobjectB <- prepare_mofa(object=MOFAobjectB, data_options=data_opts, model_options=model_opts, training_options=train_opts)
save(MOFAobjectB, file = "MOFAobjectB.RData")

## Run MOFA and save latent factors 
MOFAobject.trained <- run_mofa(MOFAobjectB, save_data = T, outfile = "/.../MOFAobject_trained.hdf5")

# Trained MOFA with the following characteristics: 
# Number of views: 4 
# Views names: RNA Meth Alt CNV 
# Number of features (per view): 5000 5000 163 8 
# Number of groups: 1 
# Groups names: group1 
# Number of samples (per group): 319 
# Number of factors: 10 

LFs <- as.data.frame(get_factors(MOFAobject.trained)$group1)
write.table(LFs, paste0("LFs_", nrow(LFs),".Cordin.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep="\t")

## Plots (standard MOFA plots)
# Data overview
pdf("Data_overview.pdf", w=5, h=4)
plot_data_overview_pretty(MOFAobjectB)
dev.off()

# Variance explained
pdf("Variance_explained.pdf", w=5, h=4)
plot_variance_explained(MOFAobject.trained)
dev.off()

get_variance_explained(MOFAobject.trained)
#         RNA        Meth         Alt         CNV 
# 51.22364275 63.15063953  0.02933155 13.39278263 

#                 RNA       Meth         Alt         CNV
# Factor1  17.0095864 34.9083853 0.004475459 4.967921348
# Factor2  13.6629123 10.4636928 0.002203085 0.016186321
# Factor3   1.4688788  9.4852865 0.003942324 3.610533756
# Factor4   5.9966510  0.3541287 0.001737472 1.437848569 
# Factor5   5.5339700  2.0055460 0.002652486 0.001582829
# Factor6   2.3031837  0.9192887 0.002361510 1.551165614
# Factor7   0.9140548  2.4325518 0.001964884 1.160429567
# Factor8   3.2178365  0.4558134 0.001677059 0.592195901 
# Factor9   2.0112523  1.5500395 0.001877455 0.002339994 
# Factor10  1.1709367  1.6503255 0.006441824 0.176142428

# Factor correlations (within model)
pdf("Factor_correlations.pdf", w=5,h=4)
plot_factor_cor(MOFAobject.trained)
dev.off()

# Factor outliers
pdf("Factor_outliers.pdf", w=5,h=4)
plot_factor(MOFAobject.trained, 1:10, dot_size = 1, dot_alpha = 0.7)
dev.off()

################################################################################

## ParetoTI set up

#                 RNA       Meth         Alt         CNV
# Factor1  17.0095864 34.9083853 0.004475459 4.967921348
# Factor2  13.6629123 10.4636928 0.002203085 0.016186321
# Factor3   1.4688788  9.4852865 0.003942324 3.610533756
# Factor4   5.9966510  0.3541287 0.001737472 1.437848569 - technical only 
# Factor5   5.5339700  2.0055460 0.002652486 0.001582829
# Factor6   2.3031837  0.9192887 0.002361510 1.551165614
# Factor7   0.9140548  2.4325518 0.001964884 1.160429567
# Factor8   3.2178365  0.4558134 0.001677059 0.592195901 - technical only
# Factor9   2.0112523  1.5500395 0.001877455 0.002339994 - technical only
# Factor10  1.1709367  1.6503255 0.006441824 0.176142428

# Exlcude LFs only associated with technical features: LFs 4, 8 and 9 
# Re-order LFs by proportion of variance in gene expression data explained 

#                 RNA       Meth         Alt         CNV
# Factor1  17.0095864 34.9083853 0.004475459 4.967921348 = 56.89037
# Factor2  13.6629123 10.4636928 0.002203085 0.016186321 = 24.14499
# Factor5   5.5339700  2.0055460 0.002652486 0.001582829 = 7.543751
# Factor6   2.3031837  0.9192887 0.002361510 1.551165614 = 4.776
# Factor3   1.4688788  9.4852865 0.003942324 3.610533756 = 14.56864
# Factor10  1.1709367  1.6503255 0.006441824 0.176142428 = 3.003846
# Factor7   0.9140548  2.4325518 0.001964884 1.160429567 = 4.509001

# Test LFs: 1+2, 1+2+5, 1+2+5+6, 1+2+5+6+3, 1+2+5+6+3+10, 1+2+5+6+3+10+7 

# Read in LFs, exclude unwanted factors and re-order columns  
LFs <- read.table("LFs_319.Cordin.txt")
LFs <- LFs[,-c(4,8,9)] # Exclude LFs only associated with technical artefacts (q-value > 0.05)

colnames(LFs) # 1: Factor1, 2: Factor2, 3: Factor3, 4: Factor5, 5: Factor6, 6: Factor7, 7: Factor10
# order required: LF1, LF2, LF5, LF6, LF3, LF10, LF7
LFs <- LFs[,c(1, 2, 4, 5, 3, 7, 6)]
save(LFs, file="Inputs_ParetoLFs.RData")

## Run ParetoTI
library(ParetoTI)
LFs <- load("Inputs_ParetoLFs.RData")

arc_ks.LFs <- lapply(2:ncol(LFs), function(k) k_fit_pch(t(LFs[,1:k]), ks = 2:6, check_installed = T, bootstrap = T, bootstrap_N = 200, maxiter = 1000, 
                                                        bootstrap_type = "s", seed = 2543, volume_ratio = "t_ratio", delta=0, conv_crit = 1e-04, order_type = "align", sample_prop = 0.75))
save(arc_ks.LFs, file = "arc_ks.LFs.RData")

arc_random.LFs <- lapply(2:ncol(LFs), function(k) randomise_fit_pch(data=t(LFs[,1:k]), arc_data = arc_ks.LFs[[k-1]], n_rand = 1000,
                                                                    bootstrap_N = 200, maxiter = 1000, type = "m", volume_ratio = "t_ratio", 
                                                                    delta=0, conv_crit = 1e-04, order_type = "align",
                                                                    sample_prop = 0.75))
save(arc_random.LFs, file = "arc_random.LFs.RData")

################################################################################

## Analyse ParetoTI fits and select which ones to examine 

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ParetoTI)

# Load data 
LFs <- read.table("LFs_319.Cordin.txt")
LFs <- LFs[,c(1,2,5,6,3,10,7)] # Exclude LFs only associated with technical features (q-value > 0.05) and re-order by % variance explained LF1, LF2, LF5, LF6, LF3, LF10, LF7
load("arc_ks.LFs.RData")
load("arc_random.LFs.RData")
load("variables_MOFA_LNET.RData")

# Get t-ratio, variance explained and total variance for each combination (summary_table)
# required to select number of LFs and values of K
summary_table <- bind_rows(lapply(2:ncol(LFs), function(i) cbind(arc_ks.LFs[[i-1]]$summary,LF=i)))
summary_table_random <- as.data.frame(bind_rows(lapply(2:ncol(LFs), function(i) cbind(arc_random.LFs[[i-1]]$obs_dist,LF=i))))

# Plot t-ratio, variance explained and total variance 
plot_tratio <- summary_table[!(is.na(summary_table$t_ratio)),]
plot_tratio$LF <- factor(plot_tratio$LF, levels=c(1:10))
pdf("t_ratio_plot.pdf", h=4,w=5)
ggplot(data=plot_tratio, aes(x=k, y=t_ratio, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")) +
  theme_classic() + labs(title="t-ratio")
dev.off()

plot_v <- summary_table
plot_v$LF <- factor(plot_v$LF, levels=c(1:10))
pdf("variance_explained_plot.pdf", h=4,w=5)
ggplot(data=plot_v, aes(x=k, y=varexpl, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")) +
  theme_classic() + labs(title="variance explained")
dev.off()

pdf("total_variance_plot.pdf", h=4,w=5)
ggplot(data=plot_v, aes(x=k, y=total_var, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")) +
  theme_classic() + labs(title="total variance")
dev.off()

p1 <- ggplot(data=plot_tratio, aes(x=k, y=t_ratio, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),name="No. LFs") +
  theme_classic() + labs(title="t-ratio", y="t-ratio") + theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(data=plot_v, aes(x=k, y=varexpl, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),name="No. LFs") +
  theme_classic() + labs(title="variance explained", y="variance explained") + theme(plot.title = element_text(hjust = 0.5))
p3 <- ggplot(data=plot_v, aes(x=k, y=total_var, group=LF)) +
  geom_line(aes(color=LF)) +
  geom_point(aes(color=LF)) + 
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),name="No. LFs") +
  theme_classic() + labs(title="total variance", y="total variance") + theme(plot.title = element_text(hjust = 0.5))
library(patchwork)
patchwork <- p1 + p2 + p3 + plot_layout(ncol=2) +  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) &  theme(plot.tag = element_text(size = 20))
pdf("ParetoTI_QC_plots.pdf", h=6.5,w=8)
patchwork
dev.off() 

# List significant fits
summary_table_random[which(summary_table_random$p_value <= 0.05),]
#    k  var_name     var_obs p_value LF
# 4  3   varexpl 0.999522916   0.001  2 - 1
# 5  3   t_ratio 0.791150880   0.001  2 - 1
# 6  3 total_var 0.002486504   0.001  2 - 1

# 64 3   varexpl 0.535616754   0.001  6 - 12
# 65 3   t_ratio 0.559410877   0.001  6 - 2
# 66 3 total_var 0.032130640   0.017  6 - 4

# 79 3   varexpl 0.482738710   0.001  7 - 13
# 80 3   t_ratio 0.534556519   0.001  7 - 3
# 81 3 total_var 0.015765398   0.005  7 - 3


# 22 4   varexpl 0.996428934   0.001  3 - 2
# 23 4   t_ratio 0.518430918   0.001  3 - 4
# 24 4 total_var 0.002930241   0.001  3 - 2

# 37 4   varexpl 0.852984651   0.001  4 - 5
# 38 4   t_ratio 0.465070543   0.001  4 - 5
# 39 4 total_var 0.035701407   0.001  4 - 5

# 52 4   varexpl 0.736447603   0.001  5 - 9
# 53 4   t_ratio 0.397384887   0.001  5 - 6
# 54 4 total_var 0.121997407   0.001  5

# 67 4   varexpl 0.658339178   0.001  6 - 10
# 68 4   t_ratio 0.311919573   0.001  6 - 7
# 69 4 total_var 0.172391940   0.001  6

# 82 4   varexpl 0.595198234   0.001  7 - 11
# 83 4   t_ratio 0.213740622   0.001  7 - 8
# 84 4 total_var 0.257298695   0.003  7


# 40 5   varexpl 0.953105982   0.001  4 - 3
# 41 5   t_ratio 0.166722011   0.001  4 - 9
# 42 5 total_var 0.180803391   0.007  4

# 70 5   varexpl 0.762930545   0.001  6 - 8
# 71 5   t_ratio 0.052756776   0.001  6
# 72 5 total_var 0.250620514   0.001  6


# 58 6   varexpl 0.927751773   0.001  5 - 4
# 59 6   t_ratio 0.072545868   0.001  5
# 60 6 total_var 0.177840327   0.001  5

# 73 6   varexpl 0.851317302   0.001  6 - 6
# 74 6   t_ratio 0.053375298   0.001  6
# 75 6 total_var 0.173268933   0.001  6

# 88 6   varexpl 0.779121354   0.001  7 - 7
# 89 6   t_ratio 0.039152588   0.001  7
# 90 6 total_var 0.458803897   0.007  7

# Three significant fits for 3 archetypes (2, 6 or 7 LFs)
# Five significant fits for 4 archetypes (3, 4, 5, 6 or 7 LFs)
# Two significant fits for 5 archetypes (4 or 6 LFs)
# Three significant fits for 6 archetypes (5, 6 or 7 LFs)

library(openxlsx)
write.xlsx(summary_table_random, file="summary_table_random.xlsx")
# Examine k=3 over two LFs (LF1+LF2), k=4 over three LFs (LF1+LF2+LF5)

# Run no boostrapping archetypes, and bootstrapped versions specifying K and LFs 
arc_ks.LFs.noboot <- lapply(2:ncol(LFs), function(k) k_fit_pch(t(LFs[,1:k]), ks = 2:6, check_installed = T, bootstrap = F, volume_ratio = "t_ratio"))

arc_3LF2.boot <- fit_pch_bootstrap(t(LFs[,1:2]), n = 200, sample_prop = 0.75, seed = 235, noc = 3, delta = 0, conv_crit = 1e-04, type = "s")
arc_4LF3.boot <- fit_pch_bootstrap(t(LFs[,1:3]), n = 200, sample_prop = 0.75, seed = 235, noc = 4, delta = 0, conv_crit = 1e-04, type = "s")

################################################################################

# Examine K3 over LFs 1 and 2

## Get archetype positions in LF space 
load("arc_3LF2.boot.RData") # Load boostrapped archetype call
arc_3LF2.boot_best <- average_pch_fits(arc_3LF2.boot) # get the best fit
arc_3_LF2 <- as.data.frame(t(arc_3LF2.boot_best$XC))
arc_3_LF2$name <- c("Arc1","Arc2","Arc3")

# Make dataframes for drawing line segments on LF maps
segments_LF12 <- data.frame(x = arc_3_LF2$Factor1, y = arc_3_LF2$Factor2, xend = arc_3_LF2$Factor1[c(2,3,1)], yend = arc_3_LF2$Factor2[c(2,3,1)])

# Plot archetype positions in LF space
# LF1 vs LF2
pdf("ParetoFit_K3LF2_LF1_LF2.pdf", h=5,w=7)
ggplot() + xlab("Latent Factor 1") + ylab("Latent Factor 2") + 
  geom_point(data = var_LNET, aes(x = Factor1, y = Factor2, color=type), size = 2.5, alpha = 0.8) +
  scale_colour_manual(name = "Tumour type", values=c("Carcinoid"="#CEDB80","Typical"="#9DB802","Atypical"="#025B0E","NET G3"="#B89f4AFF")) +
  geom_point(data = arc_3_LF2, aes(x = Factor1, y = Factor2), size = 1.5, color = "#332288") +
  geom_segment(data = segments_LF12, aes(x = x, y = y, xend = xend, yend = yend), color = "#332288", lineend = "round", linejoin = "round", size = 0.5) +
  theme_classic() + labs(title="ParetoTI fit MOFA LNET \n k=3 over two latent factors (LFs 1 & 2)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Get archetype proportions for each sample - which archetype do they 'belong' to? 
load("arc_ks.LFs.noboot.RData") # from non-boostrapped call
proportions_k3LFs2 <- as.data.frame(t(arc_ks.LFs.noboot[[1]]$pch_fits$S[[2]])) # Here we take LFs1:2 (arc_ks.LFs.noboot[[1]]) and k=3 (pch_fits$S[[2]])

# Compare the archetype positions between non-bootstrapped and bootstrapped tests
arc_ks.LFs.noboot[[1]]$pch_fits$XC[[2]] # not bootstrapped
#              [,1]      [,2]       [,3]
# Factor1 -2.953489  3.649719 0.02507211
# Factor2 -1.417870 -1.035787 2.56210276
as.data.frame(arc_3LF2.boot_best$XC) # bootstrapped 
#                V1        V2          V3
# Factor1 -2.859513  3.632725 0.008671552
# Factor2 -1.380592 -1.027027 2.548309108

# Add official archetype for each sample
proportions_k3LFs2$archetype <- sapply(rownames(proportions_k3LFs2), function(x) colnames(proportions_k3LFs2)[which.max(proportions_k3LFs2[x,])])
# smallest proportion of each archetype for a sample to be assigned that archetype: 
min(proportions_k3LFs2$V1[which(proportions_k3LFs2$archetype=="V1")]) # 0.3827974
min(proportions_k3LFs2$V2[which(proportions_k3LFs2$archetype=="V2")]) # 0.3763256
min(proportions_k3LFs2$V3[which(proportions_k3LFs2$archetype=="V3")]) # 0.4860147

plot(sort(proportions_k3LFs2$V1, decreasing=TRUE), main="Archetype_1_proportions_k3LFs2") # Arc1_proportions_k3_LFs2_min 700x600
abline(h=0.3827974,lty=2)
plot(sort(proportions_k3LFs2$V2, decreasing=TRUE), main="Archetype_2_proportions_k3LFs2") # Arc2_proportions_k3_LFs2_min 700x600
abline(h=0.3763256,lty=2)
plot(sort(proportions_k3LFs2$V3, decreasing=TRUE), main="Archetype_3_proportions_k3LFs2") # Arc3_proportions_k3_LFs2_min 700x600
abline(h=0.4860147,lty=2)

save(proportions_k3LFs2, file="proportions_K3LFs2.RData")

## Add K3 archetype 
var_LNET$archetype_k3_LF2 <- sapply(var_LNET$sample_id, function(x) proportions_k3LFs2$archetype[which(rownames(proportions_k3LFs2)==x)])

################################################################################

## Examine K4 over LFs 1, 2 and 5

## Get archetype positions in LF space 
load("arc_4LF3.boot.RData") # Load bootstrapped archetype call
arc_4LF3.boot_best <- average_pch_fits(arc_4LF3.boot) # get the best fit
arc_4_LF3 <- as.data.frame(t(arc_4LF3.boot_best$XC))
arc_4_LF3$name <- c("Arc1","Arc2","Arc3","Arc4")

# Make dataframes for drawing line segments on LF maps
segments_LF12_v2 <- data.frame(x = arc_4_LF3$Factor1, y = arc_4_LF3$Factor2, xend = arc_4_LF3$Factor1[c(2,3,1,1)], yend = arc_4_LF3$Factor2[c(2,3,1,1)])
temp <- data.frame(x=arc_4_LF3[4,1], y=arc_4_LF3[4,2], xend=arc_4_LF3$Factor1[2], yend=arc_4_LF3$Factor2[2])
segments_LF12_v2 <- rbind(segments_LF12_v2, temp)
temp <- data.frame(x=arc_4_LF3[4,1], y=arc_4_LF3[4,2], xend=arc_4_LF3$Factor1[3], yend=arc_4_LF3$Factor2[3])
segments_LF12_v2 <- rbind(segments_LF12_v2, temp)
rm(temp)

segments_LF15 <- data.frame(x = arc_4_LF3$Factor1, y = arc_4_LF3$Factor5, xend = arc_4_LF3$Factor1[c(2,3,1,1)], yend = arc_4_LF3$Factor5[c(2,3,1,1)])
temp <- data.frame(x=arc_4_LF3[4,1], y=arc_4_LF3[4,3], xend=arc_4_LF3$Factor1[2], yend=arc_4_LF3$Factor5[2])
segments_LF15 <- rbind(segments_LF15, temp)
temp <- data.frame(x=arc_4_LF3[4,1], y=arc_4_LF3[4,3], xend=arc_4_LF3$Factor1[3], yend=arc_4_LF3$Factor5[3])
segments_LF15 <- rbind(segments_LF15, temp)
rm(temp)

segments_LF25 <- data.frame(x = arc_4_LF3$Factor2, y = arc_4_LF3$Factor5, xend = arc_4_LF3$Factor2[c(2,3,1,1)], yend = arc_4_LF3$Factor5[c(2,3,1,1)])
temp <- data.frame(x=arc_4_LF3[4,2], y=arc_4_LF3[4,3], xend=arc_4_LF3$Factor2[2], yend=arc_4_LF3$Factor5[2])
segments_LF25 <- rbind(segments_LF25, temp)
temp <- data.frame(x=arc_4_LF3[4,2], y=arc_4_LF3[4,3], xend=arc_4_LF3$Factor2[3], yend=arc_4_LF3$Factor5[3])
segments_LF25 <- rbind(segments_LF25, temp)
rm(temp)

# Plot archetype positions in LF space

pdf("ParetoFit_K4LF3_LF1_LF2.pdf", h=5,w=7)
ggplot() + xlab("Latent Factor 1") + ylab("Latent Factor 2") + 
  geom_point(data = var_LNET, aes(x = Factor1, y = Factor2, color=type), size = 2.5, alpha = 0.8) +
  scale_colour_manual(name = "Tumour type", values=c("Carcinoid"="#CEDB80","Typical"="#9DB802","Atypical"="#025B0E","NET G3"="#B89f4AFF")) +
  geom_point(data = arc_4_LF3, aes(x = Factor1, y = Factor2), size = 1.5, color = "#332288") +
  geom_segment(data = segments_LF12_v2, aes(x = x, y = y, xend = xend, yend = yend), color = "#332288", lineend = "round", linejoin = "round", size = 0.5) +
  theme_classic() + labs(title="ParetoTI fit MOFA LNET \n k=4 over three latent factors (LFs 1, 2 & 5)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("ParetoFit_K4LF3_LF1_LF5.pdf", h=5,w=7)
ggplot() + xlab("Latent Factor 1") + ylab("Latent Factor 5") + 
  geom_point(data = var_LNET, aes(x = Factor1, y = Factor5, color=type), size = 2.5, alpha = 0.8) +
  scale_colour_manual(name = "Tumour type", values=c("Carcinoid"="#CEDB80","Typical"="#9DB802","Atypical"="#025B0E","NET G3"="#B89f4AFF")) +
  geom_point(data = arc_4_LF3, aes(x = Factor1, y = Factor5), size = 1.5, color = "#332288") +
  geom_segment(data = segments_LF15, aes(x = x, y = y, xend = xend, yend = yend), color = "#332288", lineend = "round", linejoin = "round", size = 0.5) +
  theme_classic() + labs(title="ParetoTI fit MOFA LNET \n k=4 over three latent factors (LFs 1, 2 & 5)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("ParetoFit_K4LF3_LF2_LF5.pdf", h=5,w=7)
ggplot() + xlab("Latent Factor 2") + ylab("Latent Factor 5") + 
  geom_point(data = var_LNET, aes(x = Factor2, y = Factor5, color=type), size = 2.5, alpha = 0.8) +
  scale_colour_manual(name = "Tumour type", values=c("Carcinoid"="#CEDB80","Typical"="#9DB802","Atypical"="#025B0E","NET G3"="#B89f4AFF")) +
  geom_point(data = arc_4_LF3, aes(x = Factor2, y = Factor5), size = 1.5, color = "#332288") +
  geom_segment(data = segments_LF25, aes(x = x, y = y, xend = xend, yend = yend), color = "#332288", lineend = "round", linejoin = "round", size = 0.5) +
  theme_classic() + labs(title="ParetoTI fit MOFA LNET \n k=4 over three latent factors (LFs 1, 2 & 5)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Get archetype proportions for each sample - which archetype do they 'belong' to? 
proportions_k4LFs3 <- as.data.frame(t(arc_ks.LFs.noboot[[2]]$pch_fits$S[[3]])) # Here we take LFs1:3 (arc_ks.LFs.noboot[[2]]) and k=4 (pch_fits$S[[3]])

# Compare the archetype positions between non-bootstrapped and bootstrapped tests
arc_ks.LFs.noboot[[2]]$pch_fits$XC[[3]] # not bootstrapped
#               [,1]       [,2]       [,3]        [,4]
# Factor1  0.08783357 -2.874241  3.5943441 -0.02270967
# Factor2 -0.67275692 -1.417196 -0.9666195  2.55201286
# Factor5 -5.11601772  1.222717  0.9240611  0.81610072
as.data.frame(arc_4LF3.boot_best$XC) # bootstrapped 
#                 V1        V2         V3          V4
# Factor1  0.1159747 -2.808502  3.5874244 -0.00962374
# Factor2 -0.6637019 -1.362303 -0.9804810  2.54491436
# Factor5 -5.0905836  1.194102  0.9198479  0.81095444

# Add official archetype for each sample
proportions_k4LFs3$archetype <- sapply(rownames(proportions_k4LFs3), function(x) colnames(proportions_k4LFs3)[which.max(proportions_k4LFs3[x,])])
# smallest proportion of each archetype for a sample to be assigned that archetype: 
min(proportions_k4LFs3$V1[which(proportions_k4LFs3$archetype=="V1")]) # 0.4350664
min(proportions_k4LFs3$V2[which(proportions_k4LFs3$archetype=="V2")]) # 0.3315186
min(proportions_k4LFs3$V3[which(proportions_k4LFs3$archetype=="V3")]) # 0.318076
min(proportions_k4LFs3$V4[which(proportions_k4LFs3$archetype=="V4")]) # 0.4535984

plot(sort(proportions_k4LFs3$V1, decreasing=TRUE), main="Archetype_1_proportions_k4LFs3") # Arc1_proportions_k4_LFs3_min 700x600
abline(h=0.4350664,lty=2)
plot(sort(proportions_k4LFs3$V2, decreasing=TRUE), main="Archetype_2_proportions_k4LFs3") # Arc2_proportions_k4_LFs3_min 700x600
abline(h=0.3315186,lty=2)
plot(sort(proportions_k4LFs3$V3, decreasing=TRUE), main="Archetype_3_proportions_k4LFs3") # Arc3_proportions_k4_LFs3_min 700x600
abline(h=0.318076,lty=2)
plot(sort(proportions_k4LFs3$V4, decreasing=TRUE), main="Archetype_4_proportions_k4LFs3") # Arc4_proportions_k4_LFs3_min 700x600
abline(h=0.4535984,lty=2)

save(proportions_k4LFs3, file="proportions_K4LFs3.RData")

## Add K4 archetype 
var_LNET$archetype_k4_LF3 <- sapply(var_LNET$sample_id, function(x) proportions_k4LFs3$archetype[which(rownames(proportions_k4LFs3)==x)])





