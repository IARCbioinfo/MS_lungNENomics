library(umap)
library(ggplot2)
library(RColorBrewer)

###### Loading of RNAseq matrix with all samples -------------------------------------
all_data_NET = read.table("All-NET-gene_count_matrix_1pass.csv", sep = ",", header = T,row.names = "gene_id")

###### Removal of genes on sex and mitochondrial chromosomes -------------------------
genespans = read.table("/data/references/Homo_sapiens/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans", stringsAsFactors = F, row.names = 1)
dim(genespans)
genes.nosex = genespans[!(genespans$V2 %in% c("chrM", "chrX", "chrY")), ]
genenames_NET = rownames(all_data_NET)
summary(genenames_NET %in% rownames(genespans))
all_data_NET = all_data_NET[which(rownames(all_data_NET) %in% rownames(genes.nosex)), ]

###### Normalization of the read counts ----------------------------------------------
#library(DESeq2)
colData_NET = as.matrix(colnames(all_data_NET))
DESeq_object_NET = DESeqDataSetFromMatrix(countData=all_data_NET,colData = colData_NET,design = ~1)
VST_NET = varianceStabilizingTransformation(DESeq_object_NET)
rm(DESeq_object_NET)
VST_NET = assay(VST_NET)
dim(VST_NET)

###### Selection of most variable genes ----------------------------------------------
rvdm_NET     = apply(VST_NET,1,var)
ntopdm_NET   = rev(which( cumsum( sort(rvdm_NET,decreasing = T) )/sum(rvdm_NET) <= 0.5 ))[1]
selectdm_NET = order(rvdm_NET, decreasing = TRUE)[seq_len(ntopdm_NET)]
VST50_NET    = t(VST_NET[selectdm_NET,])
dim(VST50_NET)

###### Run UMAP with n-neighbors = 100 -----------------------------------------------
UMAP_model_nn100 = umap(VST50_NET, n_neighbors = 100)

###### Plot umap ---------------------------------------------------------------------
distinctive_cols_NET_group <- c("#882255","#999933","#DDCC77","#117733","#CC6677","#808080","#00CCFF","#332288","#E41A1C")
Attributes_NET_ALL <- read.table("All-samples-info.csv", header=T, sep=",", row.names = 1, stringsAsFactors = F)
UMAP_TM_100 = data.frame(UMAP_model_nn100$layout)
colnames(UMAP_TM_100)=c("UMAP1","UMAP2")
UMAP_TM_100$Group = Attributes_NET_ALL[row.names(UMAP_TM_100),]$Group

umap1=ggplot(UMAP_TM_100, aes(x = UMAP1 , y = UMAP2, color = Group))+
  geom_point(size=2, alpha=0.8)+ scale_colour_manual(values=distinctive_cols_NET_group)+ #scale_color_brewer(palette="Spectral")+
  labs(x = "UMAP dimension 1" , y = "UMAP dimension 2", color = "") + theme_bw() + theme(legend.position = "right", legend.text = element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  guides(col = guide_legend( nrow = 10)) + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_discrete(labels = c("lung NETs Ca A1", "lung NETs Ca A2", "lung NETs Ca B","lung NETs supra-carcinoids","LCNEC","SCLC","pancreas NETs","small intestine NETs","rectum NETs"))

pdf("UMAP-NET.pdf",width = 20, height = 15)
grid.arrange(umap1, ncol = 1, nrow = 1)
dev.off()
