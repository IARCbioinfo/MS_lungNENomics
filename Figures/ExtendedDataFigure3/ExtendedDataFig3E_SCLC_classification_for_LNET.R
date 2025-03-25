# Author: Emilie 
# Npte: 
# Lib ---------------------------------------------------------------------

.libPaths("/home/mathiane/R/x86_64-pc-linux-gnu-library/4.1")

library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(optimx)
library(ggExtra)
library(stringr)
library(tibble)
library(lme4)
library(gridExtra)
library(tools)
library(janitor)
library(readxl)
library(scales)
library(rstatix)

# working dir -------------------------------------------------------------

setwd("/data/lungNENomics/work/MathianE/POU2F3_ASCL1_NEUROD1_YAP1")

# Palette -----------------------------------------------------------------

source("/data/lungNENomics/work/Descriptive_manuscript_data/Colour_palettes.R")


# Data --------------------------------------------------------------------
load("/data/lungNENomics/work/Descriptive_manuscript_data/RNA_sequencing_data/gene_TPM_nosex_matrix_PCA_LCNEC_SCLC_ITH_TR.RData")
load("/data/lungNENomics/work/Descriptive_manuscript_data/MOFA/MOFA.Exp.Meth.Alt.CNV_lungNENomicsCombined/ParetoTI/variables_archetypes_MOFA_LNET.RData")
load("/data/lungNENomics/work/Descriptive_manuscript_data/gencode_hg38_v33.RData")
load("/data/lungNENomics/work/Descriptive_manuscript_data/RNA_sequencing_data/vst_expr_matrix_PCA_LCNEC_SCLC_ITH_TR.RData")

dim(TPM_nosex)
head(var_LNET)
var_LNET$sample_id

## Format
ref.33.gene$gene_id[which(ref.33.gene$gene_name=="ASCL1")]# "ENSG00000139352.4"
ASCL1_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="ASCL1")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="NEUROD1")]# "ENSG00000162992.3"
NEUROD1_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="NEUROD1")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="POU2F3")]# ENSG00000137709.10"
POU2F3_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="POU2F3")]

## Gene sclc I from
# Emerging advances in defining the molecular and therapeutic landscape of small-cell lung cancer 
# Patterns of transcription factor programs and immune pathway activation define four major subtypes of SCLC with distinct therapeutic vulnerabilities
ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD274")]# ENSG00000137709.10"
CD274_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD274")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="PDCD1")]#
PDCD1_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="PDCD1")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CTLA4")]#
CTLA4_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CTLA4")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD80")]#
CD80_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD80")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD86")]#
CD86_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD86")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD38")]#
CD38_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CD38")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="TIGIT")]#
TIGIT_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="TIGIT")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="VSIR")]# synonym VISTA
VSIR_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="VSIR")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="ICOS")]# 
ICOS_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="ICOS")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="IDO1")]# 
IDO1_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="IDO1")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="LAG3")]# 
LAG3_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="LAG3")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CCL5")]# 
CCL5_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CCL5")]

ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CXCL10")]# 
CXCL10_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="CXCL10")]


tpm <- as.data.frame(t(TPM_nosex))

## Create VST TPM table ----------------
expr_gene_OI <- data.frame("sample_id" = rownames(tpm),
                           "ASCL1" = tpm[,colnames(tpm)==ASCL1_id ], 

                           "NEUROD1" = tpm[,colnames(tpm)==NEUROD1_id ],

                           "POU2F3" = tpm[,colnames(tpm)==POU2F3_id ],

                           "CXCL10" = tpm[,colnames(tpm)==CXCL10_id ],

                           "CCL5" = tpm[,colnames(tpm)==CCL5_id ],

                           "LAG3" = tpm[,colnames(tpm)==LAG3_id ],

                           "IDO1" = tpm[,colnames(tpm)==IDO1_id ],

                           "ICOS" = tpm[,colnames(tpm)==ICOS_id ],

                           "VISTA" = tpm[,colnames(tpm)==VSIR_id ],

                           "TIGIT" = tpm[,colnames(tpm)==TIGIT_id ],

                           "CD38" = tpm[,colnames(tpm)==CD38_id ],

                           "CD86" = tpm[,colnames(tpm)==CD86_id ],

                           "CD80" = tpm[,colnames(tpm)==CD80_id ],

                           "CTLA4" = tpm[,colnames(tpm)==CTLA4_id ],

                           "PDCD1" = tpm[,colnames(tpm)==PDCD1_id ],

                           "CD274" = tpm[,colnames(tpm)==CD274_id ]
                           )

dim(expr_gene_OI)
## Get archetype of the LNET cohort -------
dim(var_LNET[is.na(var_LNET$RNAseq_batch)==F & is.na(var_LNET$archetype_k4_LF3_label)==F,]  %>% dplyr::select(sample_id, archetype_k4_LF3_label))
expr_gene_OI = merge(expr_gene_OI,
      var_LNET[is.na(var_LNET$RNAseq_batch)==F & is.na(var_LNET$archetype_k4_LF3_label)==F,]  %>% dplyr::select(sample_id, archetype_k4_LF3_label),
      by="sample_id")#

## Normalization TPM ----------------------
expr_gene_OI[,2:(ncol(expr_gene_OI)-1)] = log2(expr_gene_OI[,2:(ncol(expr_gene_OI)-1)] +1)
expr_gene_OI[,2:(ncol(expr_gene_OI)-1)] = apply(expr_gene_OI[,2:(ncol(expr_gene_OI)-1)], 2, scale)

## Immune signature score -----------------------------
expr_gene_OI$I_score <- rowMeans(expr_gene_OI[,5:(ncol(expr_gene_OI)-1)] )




# Complex heatmap ---------------------------------------------------------
library(ComplexHeatmap)
## Detailed Heatmap ---------------------------------

expr_gene_OI <- expr_gene_OI[order(expr_gene_OI$archetype),]
expr_gene_OI_mat <-  t(as.matrix(expr_gene_OI[,c(2:17, 19)]))
                                                                                                      
ha_cols = HeatmapAnnotation(
  Archetype = expr_gene_OI$archetype,
  col = list(Archetype =  c(arc4)),
  annotation_legend_param = list(
    Archetype = list(direction = "horizontal", title = NULL))
    )


ht <- Heatmap(expr_gene_OI_mat,
              top_annotation  = ha_cols,  cluster_columns  = FALSE, show_column_names = FALSE , cluster_rows  = FALSE,
              heatmap_legend_param = list( direction = "horizontal" , title = NULL )#,
              
      )


draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")


## Simplify Heatmap ---------------------------------
expr_gene_OI <- expr_gene_OI[order(expr_gene_OI$archetype),]
expr_gene_OI_mat <-  t(as.matrix(expr_gene_OI[,c(2,3,4)]))

ha_cols = HeatmapAnnotation(
  Archetype = expr_gene_OI$archetype,
  col = list(Archetype =  c(arc4)),
  annotation_legend_param = list(
    Archetype = list(direction = "horizontal", title = NULL))
)


ht <- Heatmap(expr_gene_OI_mat,
              top_annotation  = ha_cols,  cluster_columns  = TRUE, show_column_names = FALSE , cluster_rows  = FALSE,
              heatmap_legend_param = list( direction = "horizontal" , title = NULL )#,
              
)

svg("Figure_hetmap_A_N_P_I.svg", width=5, height=4)

draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")

dev.off()

## Simplify Heatmap with hierarchical clustering ---------------------------------

expr_gene_OI_mat <-  t(as.matrix(expr_gene_OI[,c(2,3,4, 19)]))

ha_cols = HeatmapAnnotation(
  Archetype = expr_gene_OI$archetype,
  col = list(Archetype =  c(arc4)),
  annotation_legend_param = list(
    Archetype = list(direction = "horizontal", title = NULL))
)


ht <- Heatmap(expr_gene_OI_mat,
              top_annotation  = ha_cols,  cluster_columns  = TRUE, show_column_names = FALSE , cluster_rows  = FALSE,
              heatmap_legend_param = list( direction = "horizontal" , title = NULL )#,
              
)

draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")


# Draw box plots ----------------------------------------------------------

# ## Format data ----------------------------------------------------------


vst <- as.data.frame(t(vst_expr_matrix))
expr_gene_tmp_vst <- data.frame("sample_id" = rownames(vst),
                           
                           "ASCL1_vst" = vst[,colnames(vst) == ASCL1_id],
                           "ASCL1_tpm" = tpm[,colnames(tpm)==ASCL1_id ], 

                           "NEUROD1_vst" = vst[,colnames(vst) == NEUROD1_id],
                           "NEUROD1_tpm" = tpm[,colnames(tpm)==NEUROD1_id ],

                           "POU2F3_vst" = vst[,colnames(vst) == POU2F3_id],
                           "POU2F3_tpm" = tpm[,colnames(tpm)==POU2F3_id ],

                           "CXCL10_vst" = vst[,colnames(vst) == CXCL10_id],
                           "CXCL10_tpm" = tpm[,colnames(tpm)==CXCL10_id ],

                           "CCL5_vst" = vst[,colnames(vst) == CCL5_id],
                           "CCL5_tpm" = tpm[,colnames(tpm)==CCL5_id ],

                           "LAG3_vst" = vst[,colnames(vst) == LAG3_id],
                           "LAG3_tpm" = tpm[,colnames(tpm)==LAG3_id ],

                           "IDO1_vst" = vst[,colnames(vst) == IDO1_id],
                           "IDO1_tpm" = tpm[,colnames(tpm)==IDO1_id ],

                           "ICOS_vst" = vst[,colnames(vst) == ICOS_id],
                           "ICOS_tpm" = tpm[,colnames(tpm)==ICOS_id ],

                           "VISTA_vst" = vst[,colnames(vst) == VSIR_id],
                           "VISTA_tpm" = tpm[,colnames(tpm)==VSIR_id ],

                           "TIGIT_vst" = vst[,colnames(vst) == TIGIT_id],
                           "TIGIT_tpm" = tpm[,colnames(tpm)==TIGIT_id ],

                           "CD38_vst" = vst[,colnames(vst) == CD38_id],
                           "CD38_tpm" = tpm[,colnames(tpm)==CD38_id ],

                           "CD86_vst" = vst[,colnames(vst) == CD86_id],
                           "CD86_tpm" = tpm[,colnames(tpm)==CD86_id ],

                           "CD80_vst" = vst[,colnames(vst) == CD80_id],
                           "CD80_tpm" = tpm[,colnames(tpm)==CD80_id ],

                           
                           "CTLA4_vst" = vst[,colnames(vst) == CTLA4_id],
                           "CTLA4_tpm" = tpm[,colnames(tpm)==CTLA4_id ],

                           "PDCD1_vst" = vst[,colnames(vst) == PDCD1_id],
                           "PDCD1_tpm" = tpm[,colnames(tpm)==PDCD1_id ],

                           "CD274_vst" = vst[,colnames(vst) == CD274_id],
                           "CD274_tpm" = tpm[,colnames(tpm)==CD274_id ]

                           
                           )


dim(var_LNET[is.na(var_LNET$RNAseq_batch)==F & is.na(var_LNET$archetype_k4_LF3_label)==F,]  %>% dplyr::select(sample_id, archetype_k4_LF3_label))
expr_gene_tmp_vst = merge(expr_gene_tmp_vst,
                     var_LNET[is.na(var_LNET$RNAseq_batch)==F & is.na(var_LNET$archetype_k4_LF3_label)==F,]  %>% dplyr::select(sample_id, archetype_k4_LF3_label),
                     by="sample_id")#
# 
expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] <- log2(expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] + 1)
# expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] <- apply(expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] + 1, 2, scaled)

## Immune signature score -----------------------------
expr_gene_tmp_vst$I_score_tmp <- rowMeans(expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))])
expr_gene_tmp_vst$I_score_vst <- rowMeans(expr_gene_tmp_vst[, grep("_vst$", colnames(expr_gene_tmp_vst))])


p_ascl1_archetype <- ggplot(data=expr_gene_tmp_vst, aes(x=archetype_k4_LF3_label, y=ASCL1_tpm, fill=archetype_k4_LF3_label)) +
  geom_violin() + 
  geom_quasirandom(shape=21, col="white") +
  scale_fill_manual(values=c(arc4)) +
  stat_compare_means(
    comparisons=list(
      c("Ca A1", "Ca A2"),
      c("Ca A1", "Ca B"),
      c("Ca A1", "sc-enriched"),
      c("Ca A2", "Ca B"),
      c("Ca A2", "sc-enriched"),
      c("Ca B", "sc-enriched")
    ),
    method="t.test",
    label="p.signif",
    hide.ns=TRUE,
    tip.length=0,
    step.increase=0.07,
    
    mapping=aes(y=ASCL1_vst) # Specify ASCL1_vst for statistical comparisons
  ) +
  stat_compare_means(method="anova", label.y=max(expr_gene_tmp_vst$ASCL1_tpm, na.rm=TRUE)+12, mapping=aes(y=ASCL1_vst)) + # Use ASCL1_vst for ANOVA
  theme_classic() +
  ylab("ASCL1") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank())



p_neurod1_archetype <- ggplot(data=expr_gene_tmp_vst, aes(x=archetype_k4_LF3_label, y=NEUROD1_tpm, fill=archetype_k4_LF3_label)) +
  geom_violin() + 
  geom_quasirandom(shape=21, col="white") +
  scale_fill_manual(values=c(arc4)) +
  stat_compare_means(
    comparisons=list(
      c("Ca A1", "Ca A2"),
      c("Ca A1", "Ca B"),
      c("Ca A1", "sc-enriched"),
      c("Ca A2", "Ca B"),
      c("Ca A2", "sc-enriched"),
      c("Ca B", "sc-enriched")
    ),
    method="t.test",
    label="p.signif",
    hide.ns=TRUE,
    tip.length=0,
    step.increase=0.07,
    mapping=aes(y=NEUROD1_vst) # Specify NEUROD1_vst for statistical comparisons
  ) +
  stat_compare_means(method="anova", label.y=max(expr_gene_tmp_vst$NEUROD1_tpm, na.rm=TRUE)+12,  mapping=aes(y=NEUROD1_vst)) + # Use NEUROD1_vst for ANOVA
  theme_classic() +
  ylab("NEUROD1") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank())

p_POU2F3_archetype <- ggplot(data=expr_gene_tmp_vst, aes(x=archetype_k4_LF3_label, y=POU2F3_tpm, fill=archetype_k4_LF3_label)) +
  geom_violin() + 
  geom_quasirandom(shape=21, col="white") +
  scale_fill_manual(values=c(arc4)) +
  stat_compare_means(
    comparisons=list(
      c("Ca A1", "Ca A2"),
      c("Ca A1", "Ca B"),
      c("Ca A1", "sc-enriched"),
      c("Ca A2", "Ca B"),
      c("Ca A2", "sc-enriched"),
      c("Ca B", "sc-enriched")
    ),
    method="t.test",
    label="p.signif",
    hide.ns=TRUE,
    tip.length=0,
    step.increase=0.07,
    mapping=aes(y=POU2F3_vst) # Specify POU2F3_vst for statistical comparisons
  ) +
  stat_compare_means(method="anova", label.y=max(expr_gene_tmp_vst$POU2F3_tpm, na.rm=TRUE)+10, mapping=aes(y=POU2F3_vst)) + # Use POU2F3_vst for ANOVA
  theme_classic() +
  ylab("POU2F3") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank())



p_I_score_archetype <- ggplot(data=expr_gene_tmp_vst, aes(x=archetype_k4_LF3_label, y=I_score_tmp, fill=archetype_k4_LF3_label)) +
  geom_violin() + 
  geom_quasirandom(shape=21, col="white") +
  scale_fill_manual(values=c(arc4)) +
  stat_compare_means(
    comparisons=list(
      c("Ca A1", "Ca A2"),
      c("Ca A1", "Ca B"),
      c("Ca A1", "sc-enriched"),
      c("Ca A2", "Ca B"),
      c("Ca A2", "sc-enriched"),
      c("Ca B", "sc-enriched")
    ),
    method="t.test",
    label="p.signif",
    hide.ns=TRUE,
    tip.length=0,
    step.increase=0.07,
    mapping=aes(y=I_score_vst) # Specify I_score_vst for statistical comparisons
  ) +
  stat_compare_means(method="anova", label.y=max(expr_gene_tmp_vst$I_score_tmp, na.rm=TRUE)+10, mapping=aes(y=I_score_vst)) + # Use I_score_vst for ANOVA
  theme_classic() +
  ylab("Immune score") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank())+guides(fill=guide_legend(nrow=2))

svg("Figure_boxplot_A_N_P_I_log.svg", width=8, height=8)
ggarrange(p_ascl1_archetype, p_neurod1_archetype, p_POU2F3_archetype, p_I_score_archetype, common.legend = T, legend="bottom")
dev.off()




svg("Figure_boxplot_I_log.svg", width=3, height=5)
ggarrange( p_I_score_archetype)
dev.off()



# DLL3 --------------------------------------------------------------------


DLL3_id = ref.33.gene$gene_id[which(ref.33.gene$gene_name=="DLL3")]

vst <- as.data.frame(t(vst_expr_matrix))
expr_gene_tmp_vst <- data.frame("sample_id" = rownames(vst),
                                
                                "DLL3_vst" = vst[,colnames(vst) == DLL3_id],
                                "DLL3_tpm" = tpm[,colnames(tpm)==DLL3_id ] 
                               )
expr_gene_tmp_vst = merge(expr_gene_tmp_vst,
                          var_LNET[is.na(var_LNET$RNAseq_batch)==F & is.na(var_LNET$archetype_k4_LF3_label)==F,]  %>% dplyr::select(sample_id, archetype_k4_LF3_label),
                          by="sample_id")#
expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] <- log2(expr_gene_tmp_vst[, grep("_tpm$", colnames(expr_gene_tmp_vst))] + 1)


p_DLL3_archetype <- ggplot(data=expr_gene_tmp_vst, aes(x=archetype_k4_LF3_label, y=DLL3_tpm, fill=archetype_k4_LF3_label)) +
  geom_violin() + 
  geom_quasirandom(shape=21, col="white") +
  scale_fill_manual(values=c(arc4)) +
  stat_compare_means(
    comparisons=list(
      c("Ca A1", "Ca A2"),
      c("Ca A1", "Ca B"),
      c("Ca A1", "sc-enriched"),
      c("Ca A2", "Ca B"),
      c("Ca A2", "sc-enriched"),
      c("Ca B", "sc-enriched")
    ),
    method="t.test",
    label="p.signif",
    hide.ns=TRUE,
    tip.length=0,
    step.increase=0.07,
    mapping=aes(y=DLL3_vst) # Specify DLL3_vst for statistical comparisons
  ) +
  stat_compare_means(method="anova", label.y=max(expr_gene_tmp_vst$DLL3_tpm, na.rm=TRUE)+12, mapping=aes(y=DLL3_vst)) + # Use DLL3_vst for ANOVA
  theme_classic() +
  ylab("DLL3") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x=element_blank())+guides(fill=guide_legend(nrow=2))


svg("Figure_boxplot_DLL3_log.svg", width=3, height=5)
ggarrange( p_DLL3_archetype)
dev.off()
