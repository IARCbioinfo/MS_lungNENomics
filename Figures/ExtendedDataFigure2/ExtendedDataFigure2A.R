# Script comparison PCAWG/lungNENomics
library(tidyverse)
library(scales)
library(readxl)
library(patchwork)
library(ggrepel)
library(ggsignif)
library(ggbeeswarm)

# colors
arc4 <- c("CaA1"="#999933", "CaA2"="#DDCC77", "CaB"="#117733", "Supra-ca"="#CC6677")
pie(rep(1,4), col=arc4, labels=arc4)
type5 <- c("Typical"="#9DB802","Atypical"="#025B0E","undetermined"="#CEDB80","LCNEC"="#824833","SCLC"="#000000")
pie(rep(1, 5), col = type5, labels = type5)

## load data
### sample info
pcawg_samples = read_tsv("/data/lungNENomics/work/alcalan/Data/pcawg_sample_sheet.v1.4.2016-09-14.tsv")
pcawg_specimen = read_xlsx("/data/lungNENomics/work/alcalan/Data/pcawg_specimen_histology_August2016_v9.xlsx")

pcawg_info = left_join(pcawg_samples,pcawg_specimen) %>% filter(library_strategy=="WGS")

### find types with >=30 samples
types.tokeep = names(which(table(pcawg_info$histology_abbreviation[pcawg_info$donor_wgs_exclusion_white_gray!="Excluded"])>=30)) 


## load clinical data and archetypes
load("/data/lungNENomics/work/SextonoatesA/MOFA/Expression_Methylation/MOFA.RNA.Meth_lungNENomicsCombined/var_cli_archetypes_318.RData")
gigascience_attributes = read_tsv("/data/lungNENomics/work/SextonoatesA/Attributes.txt")

table(var_cli$archetype_k4_LF3_reordered,var_cli$molecular_group) # see correspondance with molecular groups from Alcala et al. 2019
table(var_cli$archetype_k4_LF3_reordered,var_cli$sex,var_cli$localisation_corrected) # check cluster names. B is males, A1 is F distal, A2 is young proximal

var_cli$group = var_cli$archetype_k4_LF3_reordered
var_cli$group[var_cli$archetype_k4_LF3_reordered=="V1"] = "Supra-ca"
var_cli$group[var_cli$archetype_k4_LF3_reordered=="V2"] = "CaA1"
var_cli$group[var_cli$archetype_k4_LF3_reordered=="V3"] = "CaB"
var_cli$group[var_cli$archetype_k4_LF3_reordered=="V4"] = "CaA2"

table(var_cli$group,var_cli$molecular_group)

var_cli = bind_rows(var_cli , gigascience_attributes %>% filter(!Sample_ID%in%var_cli$sample_id) %>% dplyr::rename(sample_id="Sample_ID", sex="Sex"  ) %>% mutate(age_cat = cut(as.numeric(Age),c(15.9,40.7,65.3,90.1))) )
var_cli$group[which(var_cli$Histopathology_simplified=="LCNEC")] = "LCNEC"
var_cli$group[which(var_cli$Histopathology_simplified=="SCLC")] = "SCLC"

var_cli = bind_rows(var_cli, tibble(sample_id=c("LCNEC3T","LCNEC4T","S01576","S00715"), group= "LCNEC"))

### SVs
#### PCAWG
sv.tcga.files = list.files("/data/lungNENomics/work/alcalan/Data/final_consensus_sv_bedpe_passonly.tcga.public/tcga/open/",pattern = "gz$",full.names = T)
sv.tcga.names = str_remove(list.files("/data/lungNENomics/work/alcalan/Data/final_consensus_sv_bedpe_passonly.tcga.public/tcga/open/",pattern = "gz$",full.names = F),".pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz")
sv.icgc.files = list.files("/data/lungNENomics/work/alcalan/Data/final_consensus_sv_bedpe_passonly.icgc.public/icgc/open/",pattern = "gz$",full.names = T)
sv.icgc.names = str_remove(list.files("/data/lungNENomics/work/alcalan/Data/final_consensus_sv_bedpe_passonly.icgc.public/icgc/open/",pattern = "gz$",full.names = F),".pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz")

svs = lapply(c(sv.tcga.files,sv.icgc.files) , function(x) read_tsv(gzfile(x)) )

# check if unique
table( sapply(1:length(svs), function(i) any(duplicated(svs[[i]]$sv_id)) ) ) #all unique ids
table( sapply(1:length(svs), function(i) any(duplicated(svs[[i]][,c(1:6,9:11)])) ) ) #barely any duplications

svs.num = tibble(ID=c(sv.tcga.names,sv.icgc.names), n.SVs = sapply(svs,nrow) )
svs.num = left_join(svs.num,pcawg_info,by=c("ID"="aliquot_id")) %>% filter(donor_wgs_exclusion_white_gray!="Excluded")

## lost some samples because other histology...
table(pcawg_info$histology_abbreviation) # 35 Biliary-AdenoCA

lostBACA = pcawg_info[which(pcawg_info$histology_abbreviation=="Biliary-AdenoCA"),][sapply(pcawg_info$aliquot_id[which(pcawg_info$histology_abbreviation=="Biliary-AdenoCA")], function(x) !x %in% svs.num$ID[svs.num$histology_abbreviation=="Biliary-AdenoCA"] ),] %>% 
  dplyr::select(donor_wgs_exclusion_white_gray,aliquot_id,library_strategy,project_code,histology_abbreviation,specimen_donor_treatment_type, donor_wgs_included_excluded, specimen_library_strategy)

lostBACA$aliquot_id %in% c(sv.tcga.names,sv.icgc.names) # not in files
#lostBACA$aliquot_id %in% pca_svs_by_sample2$sample # also not in files

## one sample not classified as B-ACA...
id_changed_hist = pca_svs_by_sample2$sample[which(pca_svs_by_sample2$Cancer.Types=="Biliary-AdenoCA")][!pca_svs_by_sample2$sample[which(pca_svs_by_sample2$Cancer.Types=="Biliary-AdenoCA")] %in% pcawg_info[which(pcawg_info$histology_abbreviation=="Biliary-AdenoCA"),]$aliquot_id ]
pcawg_info %>% filter(aliquot_id==id_changed_hist)

all( sapply(pca_svs_by_sample2$sample, function(x) x%in% pcawg_info$aliquot_id) ) # all there...
##

svs.num.med.per.Type = svs.num %>% group_by(histology_abbreviation) %>% summarize(med.n.SVs = median(n.SVs,na.rm=T))

#### MESOMICS
svs.num.mesomics = read_xlsx("/data/mesomics/work/mesomics1/MS/Nat_genet/R1/SI/TableS11_SVB.xlsx",skip = 2,col_types = c("text","numeric","numeric"),n_max = 115)
svs.num.mesomics.med = median(svs.num.mesomics$`Structural variants burden`)

#### lungNENomics
svs.num.lungNENomics = read_tsv("/data/lungNENomics/work/alcalan/WGS/structural_variants/sv-somatic-cns-nf_lungNENomicsAndPublic_20092022_results/SV_burden.tsv")
svs.num.lungNENomics = svs.num.lungNENomics %>% mutate(Sample_name=str_remove(Sample_name,"_T$"))

svs.num.lungNENomics = left_join(svs.num.lungNENomics,var_cli,by=c("Sample_name"="sample_id"))  %>% 
  mutate(group=factor(group,levels=c("CaA1","CaA2","CaB","Supra-ca","LCNEC")))%>% filter(!is.na(group))

svs.num.lungNENomics$consensus_pathology[svs.num.lungNENomics$group=="LCNEC"] = "LCNEC"
svs.num.lungNENomics = svs.num.lungNENomics %>% mutate(consensus_pathology =factor(consensus_pathology,levels = c("Typical","undetermined","Atypical","LCNEC")) )

svs.num.lungNENomics.med = svs.num.lungNENomics %>% group_by(group) %>% summarize(median(SVB.somatic))


library(ggbeeswarm)

summary( lm(log(SVB.somatic+1,10)~group, data=svs.num.lungNENomics )) # all different from A1 except supra-ca

gg_SVB = ggplot(  svs.num.lungNENomics , aes(x=group,y=SVB.somatic,fill=group)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=c(arc4,type5)) +
  theme_classic() + xlab("Cluster") + ylab("SVB (# SVs)") + 
  geom_signif(comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 2.6)+
  geom_signif(comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 2.8)+
  geom_signif(comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 3) + guides(fill="none")

gg_SVB.types = ggplot( svs.num.lungNENomics, aes(x=consensus_pathology,y=SVB.somatic,fill=consensus_pathology)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=type5) +
  theme_classic() + xlab("Histopathological type") + ylab("SVB (# SVs)") + 
  geom_signif(comparisons = list(c("Typical", "Atypical")), map_signif_level=TRUE,y_position = 3)+
  geom_signif(comparisons = list(c("Typical", "LCNEC")), map_signif_level=TRUE,y_position = 2.8)+
  guides(fill="none")

gg_SVB_groups_types = ggplot(  svs.num.lungNENomics , aes(x=group,y=SVB.somatic,fill=group)) + 
  geom_violin() + geom_quasirandom(aes(fill=consensus_pathology),col="white",shape=21) + scale_y_log10() + scale_fill_manual(values=c(arc4,type5)) +
  scale_color_manual(values=c(type5)) +
  theme_classic() + xlab("Cluster") + ylab("SVB (# SVs)") + 
  geom_signif(comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 2.6)+
  geom_signif(comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 2.8)+
  geom_signif(comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 3) + guides(fill="none")

# combined
summary( lm(log(SVB.somatic+1,10)~group+consensus_pathology, data=svs.num.lungNENomics %>% filter(group!="LCNEC") )) # all different from A1 except supra-ca


### SNVs
#### PCAWG
#snv = read_tsv("/data/mesomics/work/alcalan/Data/final_consensus_passonly.snv_mnv_indel.icgc.public.maf/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf")
#snvs.num = snv %>% group_by(Tumor_Sample_Barcode) %>% summarize(n.SNVs = n())
### use mut sigs instead
pcawg_sigs=read.table("/data/mesomics/files/mesomics1/CNV_Calling/2022-04-28/scripts/unsorted/Rnotebooks/mesomics/manuscript_desc/data/PCAWG/Signatures/SNV/PCAWG_sigProfiler_SBS_signatures_in_samples.csv.gz",h=T,sep=",")
pcawg_sigs$nmuts = rowSums(pcawg_sigs[,-(1:3)])
pcawg_sigs2 = left_join(pcawg_sigs,pcawg_info,by=c("Sample.Names"="icgc_specimen_id")) %>% filter(donor_wgs_exclusion_white_gray!="Excluded")
snvs.num.med.per.Type.sigs = pcawg_sigs2 %>% group_by(histology_abbreviation) %>% summarize(med.n.SNVs = median(nmuts,na.rm=T))

#### lungNENomics
snv.lungNENomics.SBS = read_tsv("/data/lungNENomics/work/alcalan/WGS/small_variants/somatic/release2_intersectmnps_15042023/mutational_signatures/input/output/SBS/input.SBS96.all")
snv.lungNENomics.SBS = tibble(Sample_name = colnames(snv.lungNENomics.SBS)[-1], n.SBS = colSums(snv.lungNENomics.SBS[,-1] ) )
snv.lungNENomics.ID  = read_tsv("/data/lungNENomics/work/alcalan/WGS/small_variants/somatic/release2_intersectmnps_15042023/mutational_signatures/input/output/ID/input.ID83.all")
snv.lungNENomics.ID  = tibble(Sample_name = colnames(snv.lungNENomics.ID)[-1], n.ID = colSums(snv.lungNENomics.ID[,-1] ) )
snv.lungNENomics = left_join(snv.lungNENomics.SBS,snv.lungNENomics.ID) %>% mutate(Sample_name = str_replace(str_remove(Sample_name,"_final"),"-","_") )

vars.lungNENomics = left_join(svs.num.lungNENomics,snv.lungNENomics) #snv.lungNENomics %>% filter(Sample_name%in% var_cli$sample_id) %>% group_by(sample_id) %>% summarise(n.SNVs=n())

snv.num.lungNENomics.med = vars.lungNENomics %>% group_by(group) %>% summarize(med.n.SV =median(SVB.somatic),med.n.SBS=median(n.SBS),med.n.ID=median(n.ID))

#svs_snvs.num.med.lungNENomics = tibble(histology_abbreviation=svs.num.lungNENomics.med$archetype_k4[1:4], 
#                                       med.n.SVs = svs.num.lungNENomics.med$`median(SVB.somatic)`[1:4] , 
#                                   med.n.SNVs = (left_join(snvs.num.lungNENomics,var_cli) %>% group_by(archetype_k4) %>% summarize(med.n.SNVs=median(n.SNVs)))[[2]] )
#svs_snvs.num.med.per.Type = bind_rows(left_join(snvs.num.med.per.Type.sigs,svs.num.med.per.Type),svs_snvs.num.med.lungNENomics)

### TMB comparison
gg_TMB = ggplot( vars.lungNENomics  , aes(x=group,y=n.SBS,fill=group)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=c(arc4,type5) ) +
  theme_classic() + xlab("Cluster") + ylab("TMB (# SNVs)") + 
  geom_signif(comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 4.3)+
  geom_signif(comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 4.4)+
  geom_signif(comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 4.5) + guides(fill="none")

gg_TMB.ID = ggplot( vars.lungNENomics  , aes(x=group,y=n.ID,fill=group)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=c(arc4,type5) ) +
  theme_classic() + xlab("Cluster") + ylab("TMB (# indels)") + 
  geom_signif(comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 2.8)+
  geom_signif(comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 2.9)+
  geom_signif(comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 3.0) + guides(fill="none")

gg_TMB.types = ggplot( vars.lungNENomics  ,aes(x=consensus_pathology,y=n.SBS,fill=consensus_pathology)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=type5) +
  theme_classic() + xlab("Histopathological type") + ylab("TMB (# SNVs)") + 
  geom_signif(comparisons = list(c("Typical", "Atypical")), map_signif_level=TRUE,y_position = 4.3)+
  guides(fill="none")

gg_TMB.ID.types = ggplot( vars.lungNENomics  ,aes(x=consensus_pathology,y=n.ID,fill=consensus_pathology)) + 
  geom_violin() + geom_quasirandom() + scale_y_log10() + scale_fill_manual(values=type5) +
  theme_classic() + xlab("Histopathological type") + ylab("TMB (# indels)") + 
  geom_signif(comparisons = list(c("Typical", "Atypical")), map_signif_level=TRUE,y_position = 2.8)+
  guides(fill="none")

ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.png",
       (gg_TMB.types+gg_TMB)/(gg_TMB.ID.types+gg_TMB.ID)/(gg_SVB.types + gg_SVB) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')), 
       width = 9,height = 4.5*3)
ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.svg",
       (gg_TMB.types+gg_TMB)/(gg_SVB.types + gg_SVB) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')),
       width = 9,height = 9)
ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.pdf",
       (gg_TMB.types+gg_TMB)/(gg_SVB.types + gg_SVB) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')), 
       width = 9,height = 9)


### CNVs
#### PCAWG
cnv.files = list.files("/data/mesomics/work/mesomics1/alcalan/PublicData/consensus.20170119.somatic.cna.annotated/",pattern = "txt$",full.names = T)
cnv.names = str_remove(list.files("/data/mesomics/work/mesomics1/alcalan/PublicData/consensus.20170119.somatic.cna.annotated/",pattern = "txt$",full.names = F),".consensus.20170119.somatic.cna.annotated.txt")

cnvs = lapply(cnv.files , function(x) read_tsv(x) )

cnvs.male = sapply(cnvs , function(x){any(x$chromosome=="Y")} )
cnvs.amp  = sapply(1:length(cnvs) , function(i){xtmp = cnvs[[i]] %>% 
  filter( (cnvs.male[i] & chromosome%in%c("X","Y") & total_cn>1) | (((!cnvs.male[i]) | !(chromosome%in%c("X","Y")) ) & (total_cn>2)) );return(sum(xtmp$end-xtmp$start,na.rm=T))} )
cnvs.del  = sapply(1:length(cnvs) , function(i){xtmp = cnvs[[i]] %>% 
  filter( (((!cnvs.male[i]) | !(chromosome%in%c("X","Y")) ) &  (total_cn<2)) |  (cnvs.male[i] & chromosome%in%c("X","Y") & total_cn<1) );return(sum(xtmp$end-xtmp$start,na.rm=T))} )

cnvs.num = tibble(ID=cnv.names, prop.amp = cnvs.amp/c(3088269832,3031042417 )[2-cnvs.male],
                  prop.del = cnvs.del/c(3088269832,3031042417 )[2-cnvs.male] )
cnvs.num = left_join(cnvs.num,pcawg_info,by=c("ID"="aliquot_id")) %>% filter(donor_wgs_exclusion_white_gray!="Excluded")

#### lungNENomics
lungNENomics.CNVs = read_tsv("/data/lungNENomics/work/alcalan/WGS/CNV/release1.1_PURPLE_19052023/lungNENomics.all.purple.cnv.somatic.19062023.tsv") %>% mutate(size=end-start+1 )
lungNENomics.CNVs.general = read_tsv("/data/lungNENomics/work/alcalan/WGS/CNV/release1.1_PURPLE_19052023/purple_summary_all.txt")

Males_samp = lungNENomics.CNVs.general$tumor_id[lungNENomics.CNVs.general$gender=="MALE"] 

CNB_cnvsegs = lungNENomics.CNVs %>% filter(!is.na(copyNumber.corrected.integer),
                                       ((!chromosome%in%c("chrX","chrY")|!sample%in%Males_samp) & 
                                          (as.numeric(copyNumber.corrected.integer)>2) ) | 
                                         (chromosome%in%c("chrX","chrY") & sample%in%Males_samp & as.numeric(copyNumber.corrected.integer)>1)) %>% 
  group_by(sample) %>% summarize(size.amp = sum(end-start+1),CNB.amp=n())
CNB_cnvsegs$prop.amp = NA
CNB_cnvsegs$prop.amp[CNB_cnvsegs$sample%in%Males_samp] = CNB_cnvsegs$size.amp[CNB_cnvsegs$sample%in%Males_samp]/3088269832
CNB_cnvsegs$prop.amp[!CNB_cnvsegs$sample%in%Males_samp] = CNB_cnvsegs$size.amp[!CNB_cnvsegs$sample%in%Males_samp]/3031042417

CNB_cnvsegs.prop.del = lungNENomics.CNVs %>% filter(!is.na(copyNumber.corrected.integer),
                                                ((!chromosome%in%c("chrX","chrY")|!sample%in%Males_samp) & 
                                                   (as.numeric(copyNumber.corrected.integer)<2) ) | 
                                                  (chromosome%in%c("chrX","chrY") & sample%in%Males_samp & as.numeric(copyNumber.corrected.integer)<1)) %>% 
  group_by(sample) %>% summarize(size.del = sum(end-start+1),CNB.del = n())

CNB_cnvsegs.prop.del$prop.del = NA
CNB_cnvsegs.prop.del$prop.del[CNB_cnvsegs.prop.del$sample%in%Males_samp]  = CNB_cnvsegs.prop.del$size.del[CNB_cnvsegs.prop.del$sample%in%Males_samp]/3088269832
CNB_cnvsegs.prop.del$prop.del[!CNB_cnvsegs.prop.del$sample%in%Males_samp] = CNB_cnvsegs.prop.del$size.del[!CNB_cnvsegs.prop.del$sample%in%Males_samp]/3031042417
CNB_cnvsegs.prop.del$sample = str_remove(CNB_cnvsegs.prop.del$sample,"_T$")

CNB_cnvsegs = full_join(CNB_cnvsegs.prop.del,CNB_cnvsegs)

cnvs.num.lungNENomics = left_join(vars.lungNENomics,CNB_cnvsegs,by=c(Sample_name="sample"))
cnvs.num.lungNENomics$prop.amp[is.na(cnvs.num.lungNENomics$prop.amp)] = 0 # NAs are actually 0s
cnvs.num.lungNENomics$prop.del[is.na(cnvs.num.lungNENomics$prop.del)] = 0 # NAs are actually 0s
cnvs.num.lungNENomics$CNB.del[is.na(cnvs.num.lungNENomics$CNB.del)] = 0 # NAs are actually 0s
cnvs.num.lungNENomics$CNB.amp[is.na(cnvs.num.lungNENomics$CNB.amp)] = 0 # NAs are actually 0s
#cnvs.num.lungNENomics = tibble(ID=cnvs.male.lungNENomics$sample, prop.amp = cnvs.amp.lungNENomics/c(3088269832,3031042417 )[2-cnvs.male.lungNENomics$male],
#                  prop.del = cnvs.del.lungNENomics/c(3088269832,3031042417 )[2-cnvs.male.lungNENomics$male] )

cnvs.prop.med.per.Type = cnvs.num.lungNENomics %>% group_by(consensus_pathology) %>% summarize(med.amp = median(prop.amp,na.rm=T),
                                                                                        n.samp.amp=sum(prop.amp!=0),
                                                                                        med.CNB.amp = median(CNB.amp),
                                                                                        med.del = median(prop.del,na.rm=T),
                                                                                        n.samp.del=sum(!is.na(prop.del)),
                                                                                        med.CNB.del = median(CNB.del))

cnvs.prop.med.per.group = cnvs.num.lungNENomics %>% group_by(group) %>% summarize(med.amp = median(prop.amp,na.rm=T),
                                                                                               n.samp.amp=sum(prop.amp!=0),
                                                                                               med.CNB.amp = median(CNB.amp),
                                                                                               med.del = median(prop.del,na.rm=T),
                                                                                               n.samp.del=sum(!is.na(prop.del)),
                                                                                               med.CNB.del = median(CNB.del))


## CNB comparisons continue here
cnvs.num.lungNENomics = cnvs.num.lungNENomics %>% mutate(WGD = cnvs.num.lungNENomics$Sample_name %in% str_replace(lungNENomics.CNVs.general$tumor_id[lungNENomics.CNVs.general$wholeGenomeDuplication],"-","_") )
cnvs.num.lungNENomics$sample_id[cnvs.num.lungNENomics$WGD] # missing "LNEN102_TU3" (ok, not right region) "LNEN105_TU2" (ok, not right region)

fisher.test(table(cnvs.num.lungNENomics$WGD,cnvs.num.lungNENomics$archetype_k4)) # nothing signif
fisher.test(table(cnvs.num.lungNENomics$WGD,cnvs.num.lungNENomics$consensus_pathology)) # nothing signif

summary( lm( cnvs.num.lungNENomics$CNB.amp~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )
summary( lm( cnvs.num.lungNENomics$CNB.del~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )
summary( lm( log(cnvs.num.lungNENomics$CNB.amp+1,10)~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )
summary( lm( log(cnvs.num.lungNENomics$CNB.del+1,10)~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )
summary( lm( cnvs.num.lungNENomics$prop.amp~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )
summary( lm( cnvs.num.lungNENomics$prop.del~cnvs.num.lungNENomics$archetype_k4*cnvs.num.lungNENomics$WGD) )

### TMB comparison
gg_CNB.amp = ggplot( cnvs.num.lungNENomics %>% filter(!is.na(archetype_k4))  , 
                 aes(x=archetype_k4,y=CNB.amp,fill=archetype_k4)) + 
  geom_violin(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,]) + 
  geom_quasirandom(col=c("black","red")[cnvs.num.lungNENomics$WGD+1]) + scale_fill_manual(values=arc4) + scale_y_log10()+
  theme_classic() + xlab("Molecular cluster") + ylab("CNB (# amplified segments)") + 
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 3.05)+
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 3.25)+
  #geom_signif(comparisons = list(c("CaA2", "CaB")), map_signif_level=TRUE,y_position = 3.05)+
  #geom_signif(comparisons = list(c("CaB", "Supra-ca")), map_signif_level=TRUE,y_position = 3.05)+
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 3.45) + guides(fill="none")

gg_CNB.amp.types = ggplot( cnvs.num.lungNENomics %>% filter(!is.na(archetype_k4))  , 
                       aes(x=consensus_pathology,y=CNB.amp,fill=consensus_pathology)) + 
  geom_violin(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,]) + 
  geom_quasirandom(col=c("black","red")[cnvs.num.lungNENomics$WGD+1]) + scale_fill_manual(values=type5) + scale_y_log10()+
  theme_classic() + xlab("Histopathological type") + ylab("CNB (# amplified segments)") + 
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("Typical", "Atypical")), map_signif_level=TRUE,y_position = 3.45)+
  guides(fill="none")

gg_CNB.del = ggplot( cnvs.num.lungNENomics %>% filter(!is.na(archetype_k4))  , 
                     aes(x=archetype_k4,y=CNB.del,fill=archetype_k4)) + 
  geom_violin() + geom_quasirandom(col=c("black","red")[cnvs.num.lungNENomics$WGD+1]) + scale_fill_manual(values=arc4) + scale_y_log10()+
  theme_classic() + xlab("Cluster") + ylab("CNB (# deleted segments)") + 
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "CaA2")), map_signif_level=TRUE,y_position = 1.6)+
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "CaB")), map_signif_level=TRUE,y_position = 1.6+0.15)+
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("CaA1", "Supra-ca")), map_signif_level=TRUE,y_position = 1.6+0.15*2) + guides(fill="none")

gg_CNB.del.types = ggplot( cnvs.num.lungNENomics %>% filter(!is.na(archetype_k4))  , 
                           aes(x=consensus_pathology,y=CNB.del,fill=consensus_pathology)) + 
  geom_violin() + geom_quasirandom(col=c("black","red")[cnvs.num.lungNENomics$WGD+1]) + scale_fill_manual(values=type5) + scale_y_log10()+
  theme_classic() + xlab("Histopathological type") + ylab("CNB (# deleted segments)") + 
  geom_signif(data = cnvs.num.lungNENomics[!cnvs.num.lungNENomics$WGD,], comparisons = list(c("Typical", "Atypical")), map_signif_level=TRUE,y_position = 1.5)+
  guides(fill="none")

## plots 
### SVs
ggplot(svs.num.med.per.Type,aes(x=histology_abbreviation,y=med.n.SVs,col=histology_abbreviation)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_classic()

### SVs vs SNVs
#### how many samples?
svs_snvs.num.med.per.Type2 = svs_snvs.num.med.per.Type %>% filter(histology_abbreviation%in%c(types.tokeep,"CaA1","CaA2","CaB","Supra-ca") )

FigS16Bleft = ggplot(svs_snvs.num.med.per.Type2,aes(x=med.n.SNVs,y=med.n.SVs,col=histology_abbreviation)) + 
  geom_point(size=(svs_snvs.num.med.per.Type2$histology_abbreviation%in%c("CaA1","CaA2","CaB","Supra-ca"))*2+1) + 
  geom_text_repel(label=svs_snvs.num.med.per.Type2$histology_abbreviation,
                  size=(svs_snvs.num.med.per.Type2$histology_abbreviation%in%c("CaA1","CaA2","CaB","Supra-ca"))*2+3,max.overlaps = 50) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c(arc4))+
  labs(title="",x="# Small variants (median)", y = "# SVs (median)")+
  theme_classic() +
  theme(legend.position = "none")

  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_log10() + ylim(c(0,250)) + theme_classic() + xlab("# SNV (median)") + ylab("# SVs (median)")

FigS16Bleft

#### CNVs 
cnvs.prop.med.per.Type2 = cnvs.prop.med.per.Type %>% filter(histology_abbreviation%in%c(types.tokeep,"CaA1","CaA2","CaB","Supra-ca") )

FigS16Bright <- ggplot(cnvs.prop.med.per.Type2 ,aes(x=med.amp*100,y=med.del*100,col=histology_abbreviation)) + 
  geom_point(size=(cnvs.prop.med.per.Type2$histology_abbreviation%in%c("CaA1","CaA2","CaB","Supra-ca"))*2+1) + 
  geom_text_repel(label=cnvs.prop.med.per.Type2$histology_abbreviation,
                  size=(cnvs.prop.med.per.Type2$histology_abbreviation%in%c("CaA1","CaA2","CaB","Supra-ca"))*2+3,max.overlaps = 10) +
  scale_color_manual(values = c(arc4))+
  labs(title="",x="% Genome amplified (median)", y = "% Genome deleted (median)")+
  theme_classic() +
  theme(legend.position = "none")

pdf("/data/lungNENomics/work/alcalan/Figures/FigS_scatter_plot_lungNENomics_pcawg_median.pdf",8,4)
FigS16Bleft + FigS16Bright
dev.off()


ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.png",
       (FigS16Bleft + FigS16Bright)/(gg_TMB.types+gg_TMB)/(gg_SVB.types + gg_SVB)/(gg_CNB.amp.types + gg_CNB.amp)/(gg_CNB.del.types + gg_CNB.del) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')), 
       width = 9,height = 18)
ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.svg",
       (FigS16Bleft + FigS16Bright)/(gg_TMB.types+gg_TMB)/(gg_SVB.types + gg_SVB)/(gg_CNB.amp.types + gg_CNB.amp)/(gg_CNB.del.types + gg_CNB.del) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')),
       width = 9,height = 18)
ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens.pdf",
       (FigS16Bleft + FigS16Bright)/(gg_TMB.types+gg_TMB)/(gg_SVB.types + gg_SVB)/(gg_CNB.amp.types + gg_CNB.amp)/(gg_CNB.del.types + gg_CNB.del) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')), 
       width = 9,height = 18)

# write tables
burden.tib = left_join(left_join(gg_TMB$data[,c("sample_id","n.SNVs")] , gg_SVB$data[,c("Sample_name","SVB.somatic")] ,by=c("sample_id"="Sample_name") ), 
          gg_CNB.amp$data[,c("sample_id","prop.del","CNB.del","prop.amp","CNB.amp")] )


write_tsv(burden.tib,file = "/data/lungNENomics/work/alcalan/Tables/Table_S_burdens.tsv")

# tests mentioned in main 
burden.wilcox = bind_rows(bind_cols(Alteration_type = "Small variants", test=c("A1 vs A2","A1 vs B","A1 vs Supra-ca","Typical vs Atypical" ),
                       as_tibble(rbind(as.numeric(wilcox.test(gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="CaA1"],gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="CaA2"])[c("statistic","p.value")]),
      as.numeric(wilcox.test(gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="CaA1"],gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="CaB"])[c("statistic","p.value")]),
      as.numeric(wilcox.test(gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="CaA1"],gg_TMB$data$n.SNVs[gg_TMB$data$archetype_k4=="Supra-ca"])[c("statistic","p.value")]),
      as.numeric(wilcox.test(gg_TMB$data$n.SNVs[gg_TMB$data$consensus_pathology=="Typical"],gg_TMB$data$n.SNVs[gg_TMB$data$consensus_pathology=="Atypical"])[c("statistic","p.value")]) ) ) ),
      bind_cols(Alteration_type = "SVs", test=c("A1 vs A2","A1 vs B","A1 vs Supra-ca","Typical vs Atypical" ),
                          as_tibble(rbind(as.numeric(wilcox.test(gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="CaA1"],gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="CaA2"])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="CaA1"],gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="CaB"])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="CaA1"],gg_SVB$data$SVB.somatic[gg_SVB$data$archetype_k4=="Supra-ca"])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_SVB$data$SVB.somatic[gg_SVB$data$consensus_pathology=="Typical"],gg_SVB$data$SVB.somatic[gg_SVB$data$consensus_pathology=="Atypical"])[c("statistic","p.value")])) ) ),

      bind_cols(Alteration_type = "CNV amplifications", test=c("A1 vs A2","A1 vs B","A1 vs Supra-ca","Typical vs Atypical" ),
                as_tibble(rbind(as.numeric(wilcox.test(gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="CaA1"&!gg_CNB.amp$data$WGD],gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="CaA2"&!gg_CNB.amp$data$WGD])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="CaA1"&!gg_CNB.amp$data$WGD],gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="CaB"&!gg_CNB.amp$data$WGD])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="CaA1"&!gg_CNB.amp$data$WGD],gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$archetype_k4=="Supra-ca"&!gg_CNB.amp$data$WGD])[c("statistic","p.value")]),
                   as.numeric(wilcox.test(gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$consensus_pathology=="Typical"&!gg_CNB.amp$data$WGD],gg_CNB.amp$data$CNB.amp[gg_CNB.amp$data$consensus_pathology=="Atypical"&!gg_CNB.amp$data$WGD])[c("statistic","p.value")]))) ),

      bind_cols(Alteration_type = "CNV deletions", test=c("A1 vs A2","A1 vs B","A1 vs Supra-ca","Typical vs Atypical" ),
                as_tibble(rbind(as.numeric(wilcox.test(gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="CaA1"&!gg_CNB.del$data$WGD],gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="CaA2"&!gg_CNB.del$data$WGD])[c("statistic","p.value")]),
                      as.numeric(wilcox.test(gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="CaA1"&!gg_CNB.del$data$WGD],gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="CaB"&!gg_CNB.del$data$WGD])[c("statistic","p.value")]),
                      as.numeric(wilcox.test(gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="CaA1"&!gg_CNB.del$data$WGD],gg_CNB.del$data$CNB.del[gg_CNB.del$data$archetype_k4=="Supra-ca"&!gg_CNB.del$data$WGD])[c("statistic","p.value")]),
                      as.numeric(wilcox.test(gg_CNB.del$data$CNB.del[gg_CNB.del$data$consensus_pathology=="Typical"&!gg_CNB.del$data$WGD],gg_CNB.del$data$CNB.del[gg_CNB.del$data$consensus_pathology=="Atypical"&!gg_CNB.del$data$WGD])[c("statistic","p.value")])) ) ) )

colnames(burden.wilcox)[3:4] = c("W.statistic","pvalue")

write_tsv(burden.tib,file = "/data/lungNENomics/work/alcalan/Tables/Table_S_burdens_tests.tsv")
# check a few examples
cnvs.num.KChRCC = cnvs.num %>% filter(histology_abbreviation=="Kidney-ChRCC",!is.na(prop.amp))
table(cnvs.num.KChRCC$prop.amp>0,cnvs.num.KChRCC$donor_wgs_exclusion_white_gray)

ids0 = cnvs.num.KChRCC %>% filter(prop.amp==0) %>% pull(aliquot_id)

for( x in cnvs[cnv.names %in% ids0 ]){
  print(table(x$total_cn))
}

## alternative : from paper SI
sample_info_S1 = read_xlsx("/data/mesomics/work/alcalan/Data/PCAWG_supplementaryTables/Supplementary Table 1.xlsx",skip=2)
table(sample_info_S1$histology_abbreviation)# 43 Kidney Ch RCC
cnvs.propdel.med.per.Type

ggplot(cnvs.prop.med.per.Type,aes(x=med.amp,y=med.del,col=histology_abbreviation)) + geom_point() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme_classic()


# Compare with Age
load("Descriptive_manuscript_data/combined_public_lungNENomics_clinical.RData")
load("Descriptive_manuscript_data/combined_public_lungNENomics_technical_data.RData")

burden.tib.clin = left_join(var_cli %>% filter(!is.na(WGS_batch)),
                            left_join(burden.tib,
                                      left_join(clinical_combined,technical_complete,by=c("manuscript_id", "lungNENomics_id")),
                                      by=c("sample_id") ),by="sample_id" ) %>% mutate(archetype_k4 = factor(archetype_k4,levels=c("CaA1","CaA2","CaB","Supra-ca")))

lm_age_smallvars = lm(n.SNVs~age_corrected + age_corrected*archetype_k4-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca") %>% mutate(archetype_k4=droplevels(archetype_k4)))
lm_age_smallvars_A1diff = lm(n.SNVs~age_corrected+age_corrected:(archetype_k4=="CaA1")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_smallvars_A2diff = lm(n.SNVs~age_corrected+age_corrected:(archetype_k4=="CaA2")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_smallvars_Bdiff  = lm(n.SNVs~age_corrected+age_corrected:(archetype_k4=="CaB")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_smallvars.ageonly = lm(n.SNVs~age_corrected-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
summary(lm_age_smallvars)
summary(lm_age_smallvars_A1diff)
summary(lm_age_smallvars_A2diff)
summary(lm_age_smallvars_Bdiff)
summary(lm_age_smallvars.ageonly)
anova(lm_age_smallvars.ageonly,lm_age_smallvars)#nodiff
anova(lm_age_smallvars_A1diff,lm_age_smallvars)#nodiff
anova(lm_age_smallvars_A2diff,lm_age_smallvars)#nodiff
anova(lm_age_smallvars_Bdiff,lm_age_smallvars) #nodiff
anova(lm_age_smallvars.ageonly,lm_age_smallvars) # nodiff => keep simplest

summary(lm_age_smallvars) 
summary(lm_age_smallvars.ageonly) # 53.5 mut/year, all groups similar
plot(lm_age_smallvars) # some of Talya's samples have higher TMB, probably due to higher seq coverage

# test glms
## poisson
glm_age_smallvars         = glm(n.SNVs~age_corrected+age_corrected:archetype_k4-1,data=burden.tib.clin, family = poisson)
glm_age_smallvars.ageonly = glm(n.SNVs~age_corrected,data=burden.tib.clin, family = poisson)
AIC(glm_age_smallvars) # better with inter
AIC(glm_age_smallvars.ageonly)
summary(glm_age_smallvars)
anova(glm_age_smallvars)

## neg binom 
burden.tib.clin4tests = burden.tib.clin %>% filter(archetype_k4!="Supra-ca")
burden.tib.clin4tests$CaA1 = (burden.tib.clin4tests$archetype_k4=="CaA1")*1
burden.tib.clin4tests$CaA2 = (burden.tib.clin4tests$archetype_k4=="CaA2")*1
burden.tib.clin4tests$CaB = (burden.tib.clin4tests$archetype_k4=="CaB")*1

library(MASS)
par(mfrow=c(1,2))
hist(burden.tib.clin$n.SNVs,xlim = c(0,20000),nclass=20)
hist(rnbinom(nrow(burden.tib.clin),size=3,mu = mean(burden.tib.clin$n.SNVs)),xlim = c(0,20000),nclass=10)
summary(mfulli <- glm.nb(n.SNVs~age_corrected+age_corrected:CaA2+age_corrected:CaB, data=burden.tib.clin4tests)) 
summary(mfull <- glm.nb(n.SNVs~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, data=burden.tib.clin4tests))
summary(mage <- glm.nb(n.SNVs~age_corrected-1, data=burden.tib.clin4tests))
summary(magei <- glm.nb(n.SNVs~age_corrected, data=burden.tib.clin4tests)) # best
summary(lmfulli <- lm(n.SNVs~age_corrected+age_corrected:CaA2+age_corrected:CaB, data=burden.tib.clin4tests))
summary(lmfull <- lm(n.SNVs~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, data=burden.tib.clin4tests))
summary(lmagei <- lm(n.SNVs~age_corrected, data=burden.tib.clin4tests))
summary(lmage  <- lm(n.SNVs~age_corrected-1, data=burden.tib.clin4tests))
AIC(mfulli)
AIC(mfull)
AIC(mage)
AIC(magei) # best
AIC(lmagei)
AIC(lmage)

cor.test( magei$y, magei$fitted.values )

glmnb.fit = predict(magei,newdata = tibble(age_corrected=seq(0,80,1)),se.fit = T)
glmnb.fit.tib = tibble(age_corrected=seq(0,80,1),n.SNVs=exp(glmnb.fit$fit))

glmnb.fit.ci.tib = tibble(age_corrected=c(seq(0,80,1),seq(80,0,-1)),
                       n.SNVs=c(exp(glmnb.fit$fit+glmnb.fit$se.fit*2),rev(exp(glmnb.fit$fit-glmnb.fit$se.fit*2) )) )

ggTMBage = ggplot(burden.tib.clin ,aes(x=age_corrected,y=n.SNVs )) + geom_point(aes(col=archetype_k4)) + scale_color_manual(values = arc4)+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,10000),expand = F)+ 
  geom_polygon(data=glmnb.fit.ci.tib,fill=rgb(0.5,0.5,0.5,0.5))+
  geom_line(data=glmnb.fit.tib)+
  theme_classic() + xlab("Age") + ylab("TMB (# small variants)")

ggTMBage
#lm 
ggplot(burden.tib.clin ,aes(x=age_corrected,y=n.SNVs)) + geom_point(aes(col=archetype_k4)) + scale_color_manual(values = arc4)+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,6000),expand = F)+ geom_smooth(method = "lm",formula = y~x-1)+
  theme_classic() + xlab("Age") + ylab("TMB (# small variants)")

# glm poisson
ggplot(burden.tib.clin ,aes(x=age_corrected,y=n.SNVs )) + geom_point(aes(col=archetype_k4)) + scale_color_manual(values = arc4)+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,10000),expand = F)+ 
  geom_smooth(method = "glm",formula = y~x, method.args = list(family = "poisson"),se = TRUE)+
  theme_classic() + xlab("Age") + ylab("TMB (# small variants)")

# lm on log
ggplot(burden.tib.clin ,aes(x=age_corrected,y=n.SNVs+1 )) + geom_point(aes(col=archetype_k4)) + scale_color_manual(values = arc4)+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(1,10000),expand = F)+ 
  geom_smooth(method = "lm",formula = y~x)+
  theme_classic() + xlab("Age") + ylab("TMB (# small variants)") + scale_y_log10()


lm_age_SVs = lm(SVB.somatic~age_corrected+age_corrected:archetype_k4-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_SVs_A1diff = lm(SVB.somatic~age_corrected+age_corrected:(archetype_k4=="CaA1")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_SVs_A2diff = lm(SVB.somatic~age_corrected+age_corrected:(archetype_k4=="CaA2")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_SVs_Bdiff  = lm(SVB.somatic~age_corrected+age_corrected:(archetype_k4=="CaB")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))
lm_age_SVs.ageonly = lm(SVB.somatic~age_corrected-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca"))

anova(lm_age_SVs.ageonly,lm_age_SVs)#nodiff
anova(lm_age_SVs.ageonly,lm_age_SVs_A1diff)#nodiff
anova(lm_age_SVs.ageonly,lm_age_SVs_A2diff)# nodiff
anova(lm_age_SVs.ageonly,lm_age_SVs_Bdiff) # diff
anova(lm_age_SVs.ageonly,lm_age_SVs.ageonly) # nodiff 
AIC(lm_age_SVs)
AIC(lm_age_SVs_A1diff)
AIC(lm_age_SVs_A2diff) # worst
AIC(lm_age_SVs_Bdiff) # best model
AIC(lm_age_SVs.ageonly)

summary(lm_age_SVs) # 0.19 mut/year Ca1 signif
summary(lm_age_SVs_A1diff) # not signif
summary(lm_age_SVs_A2diff) # not signif
summary(lm_age_SVs_Bdiff) # signif

plot(lm_age_SVs)

# NB GLM
summary(mfullip <- glm.nb(SVB.somatic~age_corrected+CaA2+CaB+age_corrected:CaA2+age_corrected:CaB, data=burden.tib.clin4tests))
summary(msomeip <- glm.nb(SVB.somatic~age_corrected+CaA2+CaB, data=burden.tib.clin4tests))
summary(mip <- glm.nb(SVB.somatic~CaA2+CaB, data=burden.tib.clin4tests)) 
summary(mfulli <- glm.nb(SVB.somatic~age_corrected+age_corrected:CaA2+age_corrected:CaB, data=burden.tib.clin4tests)) 
summary(mfull <- glm.nb(SVB.somatic~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, data=burden.tib.clin4tests))
summary(mage <- glm.nb(SVB.somatic~age_corrected-1, data=burden.tib.clin4tests))
summary(magei <- glm.nb(SVB.somatic~age_corrected, data=burden.tib.clin4tests)) # best
summary(lmfulli <- lm(SVB.somatic~age_corrected+age_corrected:CaA2+age_corrected:CaB, data=burden.tib.clin4tests))
summary(lmfull <- lm(SVB.somatic~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, data=burden.tib.clin4tests))
summary(lmagei <- lm(SVB.somatic~age_corrected, data=burden.tib.clin4tests))
summary(lmage  <- lm(SVB.somatic~age_corrected-1, data=burden.tib.clin4tests))
AIC(mfullip) 
AIC(msomeip) # best model
AIC(mip) # similar
AIC(mfulli) # similar
anova(msomeip,mip) # roughly same with age
AIC(mfull)
AIC(mage)
AIC(magei) # best
AIC(lmagei)
AIC(lmage)

cor.test( magei$y, magei$fitted.values )

newdataA1 = tibble(age_corrected=seq(min(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaA1"]),max(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaA1"]),1),
                   CaA2=0,CaB=0, archetype_k4="CaA1" )
newdataA2 = tibble(age_corrected=seq(min(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaA2"]),max(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaA2"]),1),
                   CaA2=1,CaB=0, archetype_k4="CaA2" )
newdataB  = tibble(age_corrected=seq(min(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaB"]),max(burden.tib.clin$age_corrected[burden.tib.clin$archetype_k4=="CaB"]),1),
                   CaA2=0,CaB=1, archetype_k4="CaB" )

glmnb.fit = list(predict(mip,newdata = newdataA1 ,se.fit = T) ,predict(mip,newdata = newdataA2 ,se.fit = T) , predict(mip,newdata = newdataB ,se.fit = T) ) 
glmnb.fit.tib.A1 = bind_cols(newdataA1,SVB.somatic=exp(glmnb.fit[[1]]$fit))
glmnb.fit.tib.A2 = bind_cols(newdataA2,SVB.somatic=exp(glmnb.fit[[2]]$fit))
glmnb.fit.tib.B  = bind_cols(newdataB,SVB.somatic=exp(glmnb.fit[[3]]$fit))
glmnb.fit.ci.tib.A1 = tibble(age_corrected=c(glmnb.fit.tib.A1$age_corrected,rev(glmnb.fit.tib.A1$age_corrected)),
                          SVB.somatic=c(exp(glmnb.fit[[1]]$fit+glmnb.fit[[1]]$se.fit*2),
                                        rev(exp(glmnb.fit[[1]]$fit-glmnb.fit[[1]]$se.fit*2) )),
                          CaA2 = 0,CaB = 0,archetype_k4="CaA1")
glmnb.fit.ci.tib.A2 = tibble(age_corrected=c(glmnb.fit.tib.A2$age_corrected,rev(glmnb.fit.tib.A2$age_corrected)),
                             SVB.somatic=c(exp(glmnb.fit[[2]]$fit+glmnb.fit[[2]]$se.fit*2),
                                           rev(exp(glmnb.fit[[2]]$fit-glmnb.fit[[2]]$se.fit*2) )),
                             CaA2 = 0,CaB = 0,archetype_k4="CaA2")
glmnb.fit.ci.tib.B = tibble(age_corrected=c(glmnb.fit.tib.B$age_corrected,rev(glmnb.fit.tib.B$age_corrected)),
                             SVB.somatic=c(exp(glmnb.fit[[3]]$fit+glmnb.fit[[3]]$se.fit*2),
                                           rev(exp(glmnb.fit[[3]]$fit-glmnb.fit[[3]]$se.fit*2) )),
                             CaA2 = 0,CaB = 0,archetype_k4="CaB")


ggSVBage = ggplot(burden.tib.clin ,aes(x=age_corrected,y=SVB.somatic)) + geom_point(aes(col=archetype_k4)) + 
  scale_color_manual(values = c(arc4) )+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,200),expand = F)+ 
  geom_polygon(data=glmnb.fit.ci.tib.A1,fill=rgb(0.5,0.5,0.5,0.5))+
  geom_polygon(data=glmnb.fit.ci.tib.A2,fill=rgb(0.5,0.5,0.5,0.5))+
  geom_polygon(data=glmnb.fit.ci.tib.B,fill=rgb(0.5,0.5,0.5,0.5))+
  geom_line(data=glmnb.fit.tib.A1,aes(col=archetype_k4))+
  geom_line(data=glmnb.fit.tib.A2,aes(col=archetype_k4))+
  geom_line(data=glmnb.fit.tib.B,aes(col=archetype_k4))+
  theme_classic() + xlab("Age") + ylab("SVB (# SVs)")

ggSVBage

## CNVs 
### amp 
lm_age_CNVs.amp = lm(CNB.amp~age_corrected+age_corrected:archetype_k4-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.amp_A1diff = lm(CNB.amp~age_corrected+age_corrected:(archetype_k4=="CaA1")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.amp_A2diff = lm(CNB.amp~age_corrected+age_corrected:(archetype_k4=="CaA2")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.amp_Bdiff  = lm(CNB.amp~age_corrected+age_corrected:(archetype_k4=="CaB")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.amp.ageonly = lm(CNB.amp~age_corrected-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))

anova(lm_age_CNVs.amp.ageonly,lm_age_CNVs.amp)#nodiff
anova(lm_age_CNVs.amp.ageonly,lm_age_CNVs.amp_A1diff)#diff
anova(lm_age_CNVs.amp.ageonly,lm_age_CNVs.amp_A2diff)# nodiff
anova(lm_age_CNVs.amp.ageonly,lm_age_CNVs.amp_Bdiff) # nodiff
AIC(lm_age_CNVs.amp)
AIC(lm_age_CNVs.amp_A1diff) # best model
AIC(lm_age_CNVs.amp_A2diff) # worst
AIC(lm_age_CNVs.amp_Bdiff) 
AIC(lm_age_CNVs.amp.ageonly)

summary(lm_age_CNVs.amp) # not signif
summary(lm_age_CNVs.amp_A1diff) # 8.2/year in A1 signif, only 0.7 in others
summary(lm_age_CNVs.amp_A2diff) # not signif
summary(lm_age_CNVs.amp_Bdiff) # not signif

# NB glm
summary(CNVamp.NB <- glm.nb(CNB.amp~age_corrected+age_corrected:archetype_k4-1,
                        data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca") %>% mutate(archetype_k4=droplevels(archetype_k4))))
anova(CNVamp.NB)

summary(m1 <- glm.nb(CNB.amp~age_corrected+age_corrected:CaA2+age_corrected:CaB, 
                     data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # no accumulation with age, not realistic...
summary(m1 <- glm.nb(CNB.amp~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, 
                     data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # A2 faster
summary(lm1 <- lm(CNB.amp~age_corrected+age_corrected:CaA2+age_corrected:CaB, 
                  data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # no accumulation with age, not realistic...
summary(lm1 <- lm(CNB.amp~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, 
                  data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # not signif

ggCNVampage = ggplot(burden.tib.clin%>%filter(!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]),aes(x=age_corrected,y=CNB.amp)) + geom_point(aes(col=archetype_k4)) + 
  scale_color_manual(values = c(arc4,"TRUE"= arc4[[2]], "FALSE"="gray") )+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,50),expand = F)+ 
  geom_smooth(aes(col=archetype_k4=="CaA2"),method = "lm",formula = y~x-1)+
  theme_classic() + xlab("Age") + ylab("CNB (# amplified segments)")

### del
lm_age_CNVs.del = lm(CNB.del~age_corrected+age_corrected:archetype_k4-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.del_A1diff = lm(CNB.del~age_corrected+age_corrected:(archetype_k4=="CaA1")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.del_A2diff = lm(CNB.del~age_corrected+age_corrected:(archetype_k4=="CaA2")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.del_Bdiff  = lm(CNB.del~age_corrected+age_corrected:(archetype_k4=="CaB")-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))
lm_age_CNVs.del.ageonly = lm(CNB.del~age_corrected-1,data=burden.tib.clin %>% filter(archetype_k4!="Supra-ca",sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))

anova(lm_age_CNVs.del.ageonly,lm_age_CNVs.del)#nodiff
anova(lm_age_CNVs.del.ageonly,lm_age_CNVs.del_A1diff)#nodiff
anova(lm_age_CNVs.del.ageonly,lm_age_CNVs.del_A2diff)# nodiff
anova(lm_age_CNVs.del.ageonly,lm_age_CNVs.del_Bdiff) # nodiff
AIC(lm_age_CNVs.del) # worst model
AIC(lm_age_CNVs.del_A1diff) # 
AIC(lm_age_CNVs.del_A2diff) # 
AIC(lm_age_CNVs.del_Bdiff) 
AIC(lm_age_CNVs.del.ageonly) # best model

summary(lm_age_CNVs.del) # not signif
summary(lm_age_CNVs.del_A1diff) # 8.2/year in A1 signif, only 0.7 in others
summary(lm_age_CNVs.del_A2diff) # not signif
summary(lm_age_CNVs.del_Bdiff) # not signif
summary(lm_age_CNVs.del.ageonly) #not signif

summary(m1 <- glm.nb(CNB.del~age_corrected+age_corrected:CaA2+age_corrected:CaB, 
                     data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # no accumulation with age for B
summary(m1 <- glm.nb(CNB.del~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, 
                     data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # B slower
summary(lm1 <- lm(CNB.del~age_corrected+age_corrected:CaA2+age_corrected:CaB, 
                  data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # no accumulation with age for B
summary(lm1 <- lm(CNB.del~age_corrected+age_corrected:CaA2+age_corrected:CaB-1, 
                  data=burden.tib.clin4tests %>% filter(archetype_k4!="Supra-ca",!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]))) # B slower

ggCNVdelage = ggplot(burden.tib.clin%>%filter(!sample_id%in%CNVsummary$tumor_id[CNVsummary$wholeGenomeDuplication]),aes(x=age_corrected,y=CNB.del)) + geom_point(aes(col=archetype_k4)) + 
  scale_color_manual(values = c(arc4,"TRUE"= arc4[[3]], "FALSE"="gray") )+
  coord_cartesian(xlim=c(0,max(burden.tib.clin$age_corrected,na.rm=T)),ylim=c(0,20),expand = F)+ 
  geom_smooth(aes(col=archetype_k4=="CaB"),method = "lm",formula = y~x-1)+
  theme_classic() + xlab("Age") + ylab("CNB (# deleted segments)")

ggsave("/data/lungNENomics/work/alcalan/Figures/Fig_burdens_age.png",
       (ggTMBage + ggSVBage) / (ggCNVampage + ggCNVdelage) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold')),
       height=3.5*2,width=4.5*2)
