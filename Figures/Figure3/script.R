

# Read data and format

table_S1 <- read.xlsx("Supplementary_Table_S1.xlsx", startRow = 4)
table_S1 <- table_S1[!(is.na(table_S1$`LNET.k.=.4.archetype`)), ]



# Create Figure 5a
plot_A <- as.data.frame(round(prop.table(table(table_S1$`LNET.k.=.4.archetype`, table_S1$type), margin=1), digits=2))

plot_B <- as.data.frame(round(prop.table(table(table_S1$`LNET.k.=.4.archetype`, table_S1$`sex.(omics.data.inferred)`), margin=1), digits=2))
plot_B$Var2 <- ifelse(plot_B$Var2=="F","Female","Male")

plot_C <- as.data.frame(round(prop.table(table(table_S1$`LNET.k.=.4.archetype`, table_S1$`age.(categorical)`), margin=1), digits=2))

plot_D <- as.data.frame(round(prop.table(table(table_S1$`LNET.k.=.4.archetype`, table_S1$location), margin=1), digits=2))
plot_D$Var2 <- ifelse(plot_D$Var2=="proximal","Proximal","Distal")

plot_5a <- rbind(plot_A, plot_B)
plot_5a <- rbind(plot_5a, plot_C)
plot_5a <- rbind(plot_5a, plot_D)

plot_5a$category <- c(rep("Type",16),rep("Sex",8),rep("Age",12),rep("Location",8))
plot_5a$Var1 <- factor(plot_5a$Var1, levels=c("sc-enriched","Ca B","Ca A2","Ca A1"))
plot_5a$Var2 <- factor(plot_5a$Var2, levels=c("Typical","Atypical","Carcinoid","NET G3","Female","Male","16-40","41-65","66-90","Proximal","Distal"))

ggplot(plot_5a, aes(x=Var2, y=Var1, colour=Var2, size=Freq)) +
  geom_point() +
  theme_classic() + 
  scale_color_manual(values=c("Carcinoid"="#CEDB80","Typical"="#9DB802","Atypical"="#025B0E","NET G3"="#B89F4A",
                              "Male"="#364B9A", "Female"="#A50026",
                              "16-40"="#feda8b", "41-65"="#f67e4b", "66-90"="#A50026",
                              "Proximal"="#999933", "Distal"="#AA4499")) +
  scale_size_continuous(range=c(-1,10), breaks = c(0,0.25,0.5,0.75,1)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(legend.position = "none") +
  
  geom_vline(xintercept=4.5, linetype='dashed', col = "#7f7f7f") +
  geom_vline(xintercept=6.5, linetype='dashed', col = "#7f7f7f") +
  geom_vline(xintercept=9.5, linetype='dashed', col = "#7f7f7f") +
  
  annotate(geom="text", x=2.5, y=4.4, label="Type", color="black") +
  annotate(geom="text", x=5.5, y=4.4, label="Sex", color="black") +
  annotate(geom="text", x=8, y=4.4, label="Age", color="black") +
  annotate(geom="text", x=10.5, y=4.4, label="Location", color="black") +
  
  geom_point(data=plot_5a[13,], pch=21, fill=NA, size=13*plot_5a[13,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[17,], pch=21, fill=NA, size=12*plot_5a[17,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[29,], pch=21, fill=NA, size=14*plot_5a[29,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[33,], pch=21, fill=NA, size=15*plot_5a[33,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[37,], pch=21, fill=NA, size=12*plot_5a[37,3], colour="black", stroke=1) +
  
  geom_point(data=plot_5a[14,], pch=21, fill=NA, size=13*plot_5a[14,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[18,], pch=21, fill=NA, size=12*plot_5a[18,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[30,], pch=21, fill=NA, size=14*plot_5a[30,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[26,], pch=21, fill=NA, size=15*plot_5a[26,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[42,], pch=21, fill=NA, size=11*plot_5a[42,3], colour="black", stroke=1) +
  
  geom_point(data=plot_5a[23,], pch=21, fill=NA, size=12*plot_5a[23,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[31,], pch=21, fill=NA, size=14*plot_5a[31,3], colour="black", stroke=1) +
  geom_point(data=plot_5a[35,], pch=21, fill=NA, size=14*plot_5a[35,3], colour="black", stroke=1) +
  
  labs(x="", y="Molecular group")

# Create Figure 5b
table_S1$smoking_2 <- ifelse(table_S1$smoking %in% c("current", "former"), "ever", table_S1$smoking)

smoking_data <- table_S1[which(table_S1$type %in% c("Typical","Atypical")),c(1,13,46,60)]
smoking_data <- smoking_data[!(is.na(smoking_data$smoking_2)),]
smoking_data$type <- ifelse(smoking_data$type=="Typical","grade-1","grade-2")
smoking_data$group_merge <- paste0(smoking_data$`LNET.k.=.4.archetype`, "_", smoking_data$type) 

smoking_plot_t <- as.data.frame(round(prop.table(table(smoking_data$`LNET.k.=.4.archetype`[which(smoking_data$type=="grade-1")], smoking_data$smoking_2[which(smoking_data$type=="grade-1")]), margin=1), digits = 2))
smoking_plot_a <- as.data.frame(round(prop.table(table(smoking_data$`LNET.k.=.4.archetype`[which(smoking_data$type=="grade-2")], smoking_data$smoking_2[which(smoking_data$type=="grade-2")]), margin=1), digits = 2))
smoking_plot_t$type <- "grade-1"
smoking_plot_a$type <- "grade-2"
smoking_plot_f <- rbind(smoking_plot_t, smoking_plot_a)

smoking_plot_f$type <- factor(smoking_plot_f$type, levels=c("grade-2","grade-1"))
p5b_1 <- ggplot(smoking_plot_f, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat="identity", position="fill") + 
  scale_fill_manual(values = c("never"="#364B9A", "ever"="#A50026")) +
  theme_classic() + guides(fill=guide_legend(title="Smoking status")) +
  labs(x="", y="Proportion") + 
  facet_wrap(~type) +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) + scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1)) +
  geom_segment(aes(x=2,xend=3,y=1.05,yend=1.05)) +  geom_text(aes(label = "**", y = 1.06, x=2.5)) +
  theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11,angle=45, hjust=1), axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) + 
  theme(legend.text=element_text(size=11), legend.title=element_text(size=11))


table <- read.table("Table_dnds_smoking.tsv",sep="\t", header=TRUE)
arc4 <- c("Ca A1"="#999933", "Ca A2"="#DDCC77", "Ca B"="#117733", "sc-enriched"="#CC6677") 
type5 <- c("Typical"="#9DB802","Atypical"="#025B0E","Carcinoid"="#CEDB80","LCNEC"="#824833","SCLC"="#000000")

table$group <- factor(table$group, levels=c("All", "Typical", "Atypical", "Ca A1", "Ca A2", "Ca B", "sc-enriched"))
table_v2 <- table[!(table$group %in% c("All")), ]

p5b_2 <- ggplot(table_v2, aes(x=group,y=mle,col=group)) + geom_point(show.legend = F, size=2) + 
  theme_classic() + 
  facet_grid(.~Smoking_status) + geom_segment(aes(y=cilow,yend=cihigh,x=group,xend=group),show.legend = F, size=0.65) + scale_color_manual(values = c(arc4,type5)) +
  geom_hline(yintercept=1,linetype="dashed") + ylab("Selection (dn/ds)") + xlab("") +  theme(axis.text.y=element_text(size=11), axis.text.x=element_text(size=11,angle=45, hjust=1), axis.title.x=element_text(size=11), axis.title.y=element_text(size=11)) + 
  theme(legend.text=element_text(size=11), legend.title=element_text(size=11))

patchwork <- p5b_1 + p5b_2
patchwork

















