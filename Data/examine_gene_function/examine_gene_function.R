
## Function to compare gene expression between groups and output a violin plot

# Function is called examine_gene
# Takes two inputs: a gene name (any) and a sample set (one of LNET or LNEN)
examine_gene("SEZ6", "LNEN") # example



examine_gene <- function(gene, samples) {
  
  library(ggplot2)
  library(ggbeeswarm)
  library(tidyverse)
  library(patchwork)
  
  load("gene_TPM_nosex_matrix_PCA_LCNEC_SCLC_ITH_TR.RData")
  load("vst_expr_matrix_PCA_LCNEC_SCLC_ITH_TR.RData")
  load("gencode_hg38_v33.RData") 
  
  metadata <- read.csv("Metadata.csv", skip=3)
  LNET <- metadata$sample.id[which(grepl("Expr", metadata$omics.group) & metadata$type %in% c("Atypical", "Typical", "Carcinoid", "NET G3"))]
  LNEC <- metadata$sample.id[which(metadata$type %in% c("LCNEC","SCLC"))]
  
  gene <- gene  
  gene_names <- ref.33.gene$gene_id[which(ref.33.gene$gene_name %in% gene)]
  print(paste0(gene_names," is the gene id for ", gene))
  
  group_colours <- c("Ca A1"="#999933", "Ca A2"="#DDCC77", "Ca B"="#117733", "sc-enriched"="#CC6677", "LCNEC"="#882255", "SCLC"="#332288")
  type_colours <- c("Typical"="#9DB802","Atypical"="#025B0E","LCNEC"="#882255", "SCLC"="#332288")
  
  samples <- samples 
  if(samples == "LNET"){
    
    vst_data <- as.data.frame(t(vst_expr_matrix[which(rownames(vst_expr_matrix) %in% gene_names), which(colnames(vst_expr_matrix) %in% LNET)]))
    colnames(vst_data) <- "vst"
    
    tpm_data <- as.data.frame(t(TPM_nosex[which(rownames(TPM_nosex) %in% gene_names), which(colnames(vst_expr_matrix) %in% LNET)]))
    colnames(tpm_data) <- "tpm"
    
    expr_data <- merge(vst_data, tpm_data, by="row.names")
    colnames(expr_data)[1] <- "sample_id"
    
    expr_data$group <- sapply(expr_data$sample_id, function(x) metadata$LNET.k...4.archetype[which(metadata$sample.id==x)])
    expr_data$type <- sapply(expr_data$sample_id, function(x) metadata$type[which(metadata$sample.id==x)])
    
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Typical", "Atypical"))] ~ 
                   expr_data$type[which(expr_data$type %in% c("Typical", "Atypical"))], alternative = "two.sided"))
    
    expr_data$type <- factor(expr_data$type, levels=c("Typical","Atypical"))
    
    p1 <- ggplot(data=expr_data[which(expr_data$type %in% c("Typical","Atypical")),], aes(x=type, y = log10(tpm+1), fill=type)) +
      geom_violin(scale = "width") +
      geom_boxplot(fill="white", width=0.1, outlier.size = 0.4) +
      scale_fill_manual(name = "Tumour type", values=c(type_colours)) +
      theme_classic() + ylab(paste0(gene," expression [log10(TPM+1)]")) + xlab("Tumour type") +
      theme(legend.position = "bottom") +
      guides(color=guide_legend(nrow=1,byrow=TRUE)) +
      scale_y_continuous(limits=c(0,max(log10(expr_data$tpm+1)))) +
      theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11),
            axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) +
      theme(legend.text=element_text(size=11), legend.title=element_text(size=11)) +
      guides(fill="none")

    print(paste0("lung NET molecular group ANOVA P value = ", summary(aov(vst ~ group, data = expr_data))[[1]][1,5]))
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "Ca A2"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "Ca A2"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "Ca B"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "Ca B"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "sc-enriched"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "Ca B"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "Ca B"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "sc-enriched"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca B", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca B", "sc-enriched"))], alternative = "two.sided"))
    
    p2 <- ggplot(data=expr_data, aes(x=group, y = log10(tpm+1), fill=group)) +
      geom_violin(scale = "width") + 
      geom_boxplot(fill="white", width=0.1, outlier.size = 0.4) +
      geom_quasirandom(shape=21, col="white",data=expr_data %>% filter(group =="sc-enriched")) +
      scale_fill_manual(name = "Molecular group", values=c(group_colours)) +
      theme_classic() + ylab(paste0(gene," expression [log10(TPM+1)]")) + xlab("Molecular group") +
      theme(legend.position = "bottom") +
      scale_y_continuous(limits=c(0,max(log10(expr_data$tpm+1)))) +
      guides(color=guide_legend(nrow=1,byrow=TRUE)) + 
      theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11), 
            axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) + 
      theme(legend.text=element_text(size=11), legend.title=element_text(size=11)) +
      guides(fill="none")
    
    patchwork <- p1 + p2 + plot_layout(widths = c(0.5,1))
    patchwork

  }
  else{
    
    vst_data <- as.data.frame(t(vst_expr_matrix[which(rownames(vst_expr_matrix) %in% gene_names), which(colnames(vst_expr_matrix) %in% c(LNET,LNEC))]))
    colnames(vst_data) <- "vst"
    
    tpm_data <- as.data.frame(t(TPM_nosex[which(rownames(TPM_nosex) %in% gene_names), which(colnames(vst_expr_matrix) %in% c(LNET,LNEC))]))
    colnames(tpm_data) <- "tpm"
    
    expr_data <- merge(vst_data, tpm_data, by="row.names")
    colnames(expr_data)[1] <- "sample_id"
    
    expr_data$group <- sapply(expr_data$sample_id, function(x) metadata$LNET.k...4.archetype[which(metadata$sample.id==x)])
    expr_data$type <- sapply(expr_data$sample_id, function(x) metadata$type[which(metadata$sample.id==x)])
    
    expr_data$type <- sapply(expr_data$sample_id, function(x) metadata$type[which(metadata$sample.id==x)])
    expr_data$group <- sapply(expr_data$sample_id, function(x) metadata$LNET.k...4.archetype[which(metadata$sample.id==x)])
    expr_data$group <- ifelse(is.na(expr_data$group), expr_data$type, expr_data$group)
    
    
    print(paste0("lung NEN type ANOVA P value = ", summary(aov(vst ~ type, data = expr_data[which(expr_data$type %in% c("Typical","Atypical","LCNEC","SCLC")),]))[[1]][1,5]))
    
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Typical", "Atypical"))] ~ 
             expr_data$type[which(expr_data$type %in% c("Typical", "Atypical"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Typical", "LCNEC"))] ~ 
             expr_data$type[which(expr_data$type %in% c("Typical", "LCNEC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Typical", "SCLC"))] ~ 
             expr_data$type[which(expr_data$type %in% c("Typical", "SCLC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Atypical", "LCNEC"))] ~ 
             expr_data$type[which(expr_data$type %in% c("Atypical", "LCNEC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$type %in% c("Atypical", "SCLC"))] ~ 
             expr_data$type[which(expr_data$type %in% c("Atypical", "SCLC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$type %in% c("LCNEC", "SCLC"))] ~ 
             expr_data$type[which(expr_data$type %in% c("LCNEC", "SCLC"))], alternative = "two.sided"))
    
    expr_data$type <- factor(expr_data$type, levels=c("Typical","Atypical","LCNEC","SCLC","Carcinoid","NET G3"))
    
    p1 <- ggplot(data=expr_data[which(expr_data$type %in% c("Typical","Atypical","LCNEC","SCLC")),], aes(x=type, y = log10(tpm+1), fill=type)) +
      geom_violin(scale = "width") + 
      geom_boxplot(fill="white", width=0.1, outlier.size = 0.4) +
      scale_fill_manual(name = "Tumour type", values=c(type_colours)) +
      theme_classic() + ylab(paste0(gene," expression [log10(TPM+1)]")) + xlab("Tumour type") +
      theme(legend.position = "bottom") +
      guides(color=guide_legend(nrow=1,byrow=TRUE)) + 
      scale_y_continuous(limits=c(0,max(log10(expr_data$tpm+1)))) +
      theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11), 
            axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) + 
      theme(legend.text=element_text(size=11), legend.title=element_text(size=11)) +
      guides(fill="none")
    
    
    print(paste0("lung NEN molecular group/type ANOVA P value = ", summary(aov(vst ~ group, data = expr_data))[[1]][1,5])) 
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "Ca A2"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "Ca A2"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "Ca B"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "Ca B"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "sc-enriched"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "LCNEC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "LCNEC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A1", "SCLC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A1", "SCLC"))], alternative = "two.sided"))
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "Ca B"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "Ca B"))], alternative = "two.sided")) 
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "sc-enriched"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "LCNEC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "LCNEC"))], alternative = "two.sided")) 
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca A2", "SCLC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca A2", "SCLC"))], alternative = "two.sided"))
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca B", "sc-enriched"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca B", "sc-enriched"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca B", "LCNEC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca B", "LCNEC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("Ca B", "SCLC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("Ca B", "SCLC"))], alternative = "two.sided"))
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("sc-enriched", "LCNEC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("sc-enriched", "LCNEC"))], alternative = "two.sided"))
    print(t.test(expr_data$vst[which(expr_data$group %in% c("sc-enriched", "SCLC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("sc-enriched", "SCLC"))], alternative = "two.sided"))
    
    print(t.test(expr_data$vst[which(expr_data$group %in% c("LCNEC", "SCLC"))] ~ 
             expr_data$group[which(expr_data$group %in% c("LCNEC", "SCLC"))], alternative = "two.sided"))
    
    expr_data$group <- factor(expr_data$group, levels=c("Ca A1","Ca A2","Ca B","sc-enriched","LCNEC","SCLC"))
    
    p2 <- ggplot(data=expr_data, aes(x=group, y = log10(tpm+1), fill=group)) +
      geom_violin(scale = "width") + 
      geom_boxplot(fill="white", width=0.1, outlier.size = 0.4) +
      geom_quasirandom(shape=21, col="white",data=expr_data %>% filter(group =="sc-enriched")) +
      scale_fill_manual(name = "Molecular group", values=c(group_colours)) +
      theme_classic() + ylab(paste0(gene, " expression [log10(TPM+1)]")) + xlab("Molecular group") +
      theme(legend.position = "bottom") +
      guides(color=guide_legend(nrow=1,byrow=TRUE)) + 
      scale_y_continuous(limits=c(0,max(log10(expr_data$tpm+1)))) +
      theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11), 
            axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) + 
      theme(legend.text=element_text(size=11), legend.title=element_text(size=11)) +
      guides(fill="none")
    
    patchwork <- p1 + p2 + plot_layout(widths = c(0.4,0.6))
    patchwork
    
  }
}



