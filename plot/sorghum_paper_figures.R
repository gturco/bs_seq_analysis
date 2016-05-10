library(grid)
library(ggplot2)

,


x = "/Users/gturco/Documents/Data/gene_heatmap/all.vascular.CpG.1.filtered.genes.txt"
vas1 = read.table(x, sep="\t", col.names= c("start","meth","gene"))
colnames(vas1) = c("start","meth","gene")
p <- ggplot(vas1,aes(gene,start)) + geom_tile(aes(fill=meth))
p


#### Figure 1 DNA methylation patterns across the gene body.
### 4 figues for figure

library(ggplot2)
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/shoot/all.shoot.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHH","","CHG"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation") 

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/shoot.pdf", width=5.60, height=3.2, dpi=600, units="cm")  

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/root/all.root.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHH","","CHG"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation")

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/root.pdf", width=5.60, height=3.2, dpi=600, units="cm")  

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/nonvas/all.nonvascular.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHH","","CHG"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation")

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/non.pdf", width=5.60, height=3.2, dpi=600, units="cm")  

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/vas/all.vascular.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHH","","CHG"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation")

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/vas.pdf", width=5.60, height=3.2, dpi=600, units="cm")  

### ledgend 


ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(text=element_text(size=7), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.3, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/tissue_leg.pdf", width=4.20, height=2.4, dpi=600, units="cm")  


#### Figure 2 DNA methylation patterns across tissue types.

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/all/CG_tissue.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/CG.pdf", width=5.60, height=3.2, dpi=600, units="cm")  

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/all/CHG_tissue.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/all/CHH_tissue.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = 0.3) +
  theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/CHH.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

### ledgend 

x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
levels(x$mtype)[levels(x$mtype)=="root"] <- "Root"
levels(x$mtype)[levels(x$mtype)=="shoot"] <- "Shoot"
levels(x$mtype)[levels(x$mtype)=="vascular"] <- "Vascular"
levels(x$mtype)[levels(x$mtype)=="nonvascular"] <- "Nonvascular"
ggplot(data=x, aes(x = pos, y = freq* 100 ,group=mtype, colour = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=7), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.3, "cm")) + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/new/tissue/all_leg.pdf", width=5.60, height=3.2, dpi=600, units="cm") 


### figure 3 Vascular methylation patterns in genic and intergenic regions.
### group root, shoot 

file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/all_shoot.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("CG", "CHG","CHH"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ctype, fill = ctype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_brewer(type = "div", palette = 5) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation") + ylim(0,80)


ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/shoot.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/all_root.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("CG", "CHG","CHH"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ctype, fill = ctype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_brewer(type = "div", palette = 5) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/root.pdf", width=4.20, height=2.4, dpi=600, units="cm") 



file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/all_non.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("CG", "CHG","CHH"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ctype, fill = ctype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_brewer(type = "div", palette = 5) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("")  + ylim(0,80)


ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/non.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/vas_all.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("CG", "CHG","CHH"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ctype, fill = ctype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_brewer(type = "div", palette = 5) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/vas.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

### ledgend 

levels(x$ctype)[levels(x$ctype)=="CDs"] <- "CDS"
levels(x$ctype)[levels(x$ctype)=="5 promoter"] <- "1Kb 5'"
levels(x$ctype)[levels(x$ctype)=="3 promoter"] <- "1Kb 3'"
levels(x$ctype)[levels(x$ctype)=="miRNAs"] <- "microRNAs"
levels(x$ctype)[levels(x$ctype)=="tandem duplicates"] <- "tandem "
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ctype, fill = ctype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_brewer(type = "div", palette = 5) +
  theme_classic() + theme(text=element_text(size=7), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.2, "cm")) 

25499034
### figure 3 Vascular methylation patterns in genic and intergenic regions.
###GROUP ROOT SHOOT VAS NONVAS INTO ONE FILE 

## CG
file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/CG.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
x$mtype <- factor(x$mtype, levels = c("Root","Shoot","Vascular","Nonvascular"))
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = ctype, y = freq* 100 ,group=mtype, fill = mtype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436")) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/CG.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


# CHG
file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/CHG.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
TSS
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = ctype, y = freq* 100 ,group=mtype, fill = mtype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436")) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/CHG.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


#CHH

file_name = "/Users/gturco/Documents/Data/Sorg/meth_bar/all_CHH.txt"
x = read.table(file_name, col.names=c("freq","ctype","mtype"), sep="\t", header= FALSE)
TSS
x$ctype <- factor(x$ctype, levels = c("CDs", "3 promoter","5 promoter","miRNA", "TE","tandem duplicates")) 
ggplot(data=x, aes(x = ctype, y = freq* 100 ,group=mtype, fill = mtype)) +  geom_bar(stat ="identity", position = "dodge") + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436")) + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + ylim(0,80)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_bar/CHH.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 



### Figure 4 RPKM CG

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CG.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CG.pdf",  width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.CG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/non_CG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CpG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vas_CG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme( text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm"), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.2, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#F97976","#F82F2B","#B92320","#783A39"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vas_big_CG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

### CHG

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/new/new_CHG/rpkm/all.root.filtered.CHG.txt"384
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.filtered.CHG.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/non_CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHG.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))
ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vas_CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm"), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.2, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#b9d2e2","#6aacd4","#4e91a5","#1d66a3"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vas_big_CHG.pdf", width=5.60, height=3.2, dpi=600, units="cm") 


#### CHH

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.shoot.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/shoot_CHH.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.root.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/root_CHH.pdf",width=5.60, height=3.2, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.nonvascular.filtered.CHH.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/non_CHH.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/all.vascular.CHH.filtered.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("1","2","3","4"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026"))

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/vas_CHH.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype))  +   geom_line(size = .3) + scale_color_manual(values = c("#53CA77","#15C048","#0E8130","#1A4026")) + 
  theme_classic() + theme(text=element_text(size=20))    +  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm"), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.2, "cm")) + xlab(NULL) + ylab("% methylation") 

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/leg.pdf", width=5.60, height=3.2, dpi=600, units="cm") 

#### Figure 5 RPKM over tissue types

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_1.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(14,84)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_1.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_2.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(14,84)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_2.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_3.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))

ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(14,84)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_3.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_all_4.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(14,84)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CG_4.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

####CHG
file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_all_1.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,35)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_1.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_all_2.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("")+ scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,35)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_2.pdf",width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_all_3.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,35)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_3.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_all_4.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,35)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHG_4.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


###CHH

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_all_1.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,21)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_1.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_all_2.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,21)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_2.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_all_3.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,21)

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_3.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_all_4.txt"
x = read.table(file_name, col.names=c("pos","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("root","shoot","vascular","nonvascular"))
ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("") + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8")) + ylim(5,21)
ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/CHH_4.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


### ledgend 




ggplot(data=x, aes(x = pos, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic()    +  theme(legend.key.height  = unit(.15, "cm"),text=element_text(size=7), legend.key.width  = unit(.2, "cm")) + scale_colour_manual(values=c("#c6141c","#fccb0c","#009987","#3261a8"))
ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/leg.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

ggsave("/Users/gturco/Documents/Data/Sorg/TSS/rpkm_new/tissue/Tissue/leg.pdf", width=5.0, height=2.4, dpi=600, units="cm") 



#### Figure 6 changes in methylation between tissue types like tomatopaper figure 5

file_name = "/Users/gturco/Documents/Data/Sorg/figure_5/root_non.txt"
x = read.table(file_name, col.names=c("change","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHG","CHH"))
ggplot(data=x, aes(x = change, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic()   +  theme(legend.position ="none", text=element_text(size=6), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("")

ggsave("/Users/gturco/Documents/Data/Sorg/figure_5/root_non.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/figure_5/vas_non.txt"
x = read.table(file_name, col.names=c("change","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHG","CHH"))
ggplot(data=x, aes(x = change, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=8), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("")

ggsave("/Users/gturco/Documents/Data/Sorg/figure_5/vas_non.pdf", width=4.2, height=2.4, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/figure_5/root_shoot.txt"
x = read.table(file_name, col.names=c("change","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHG","CHH"))
ggplot(data=x, aes(x = change, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=6), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% methylation")

ggsave("/Users/gturco/Documents/Data/Sorg/figure_5/root_shoot.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/Data/Sorg/figure_5/root_vas.txt"
x = read.table(file_name, col.names=c("change","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHG","CHH"))
ggplot(data=x, aes(x = change, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme(text=element_text(size=20))    +  theme(legend.position ="none", text=element_text(size=8), panel.margin = unit(0, "cm"), plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("% Methylation")

ggsave("/Users/gturco/Documents/Data/Sorg/figure_5/root_vas.pdf", width=4.2, height=2.4, dpi=600, units="cm") 

### ledgend 
ggplot(data=x, aes(x = change, y = freq * 100 ,group=mtype, color = mtype)) +  geom_line(size = .3) +
  theme_classic() + theme()  + theme(text=element_text(size=8), legend.key.height  = unit(.15, "cm"), legend.key.width  = unit(.2, "cm")) 

ggsave("/Users/gturco/Documents/Data/Sorg/figure_5/leg.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

### Average methylation using  Weighted methylation level 

file_name = "/Users/gturco/Documents/Data/Sorg/meth_avg.txt"
x = read.csv(file_name, col.names=c("ttype","freq","mtype"))
x$mtype <- factor(x$mtype, levels = c("CG","CHG","CHH"))
x$ttype <- factor(x$ttype, levels = c("Root","Shoot","Vascular","Nonvascular"))
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ttype, fill = ttype, linetype=ttype)) +  geom_bar(stat ="identity", position = "dodge") + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("Weighted Mean") + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436")) +ylim(0,60)

ggsave("/Users/gturco/Documents/Data/Sorg/meth_wavg.pdf", width=4.20, height=2.4, dpi=600, units="cm") 

### ledgend 
ggplot(data=x, aes(x = mtype, y = freq* 100 ,group=ttype, fill = ttype, linetype=ttype)) +  geom_bar(stat ="identity", position = "dodge") +
  theme_classic() + theme(text=element_text(size=20))  + theme(text=element_text(size=7), legend.key.height  = unit(.2, "cm"), legend.key.width  = unit(.2, "cm"))  + scale_fill_manual(values=c("#e1131e","#3c4a96","#009987","#71b436"))

ggsave("/Users/gturco/Documents/Data/Sorg/barleg.pdf", width=4.20, height=2.4, dpi=600, units="cm") 



coexpression

file_name = "/Users/gturco/Documents/code/Brady/sorghum/co-express/slope/coexpress_bmr2.graph_all"
x = read.table(file_name, sep="\t", header=TRUE)
levels(x$Anno)
ggplot(data=x, aes(x = tissue, y = value, group=gene, color=Anno, alpha= 1)) + geom_line() + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.05), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("RPKM Mean") + scale_colour_manual(values=c("black","#d86134", "#007742","#626362","#8d2b53","#27295d"))
ggsave("/Users/gturco/Documents/Data/Sorg/bmr2_coexprss.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


file_name = "/Users/gturco/Documents/code/Brady/sorghum/co-express/slope/coexpress_vnd7.graph_all"
x = read.table(file_name, sep="\t", header=TRUE)
levels(x$Anno)
ggplot(data=x, aes(x = tissue, y = value, group=gene, color=Anno, alpha= 1)) + geom_line() + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.05), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("RPKM Mean") + scale_colour_manual(values=c("#81213c","black","#d86134", "#777776","#27295d"))
ggsave("/Users/gturco/Documents/Data/Sorg/vnd7_coexprss.pdf", width=5.60, height=3.21, dpi=600, units="cm") 

file_name = "/Users/gturco/Documents/code/Brady/sorghum/co-express/slope/coexpress_cesa.graph_all"
x = read.table(file_name, sep="\t", header=TRUE)
levels(x$Anno)
ggplot(data=x, aes(x = tissue, y = value, group=gene, color=Anno, alpha= 1)) + geom_line() + theme_classic() + theme(legend.position ="none", text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.05), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + xlab(NULL) + ylab("RPKM Mean") + scale_colour_manual(values=c("#81213c","#872328","black","#00492d", "#777776","#27295d"))
ggsave("/Users/gturco/Documents/Data/Sorg/cesa_coexprss.pdf", width=5.60, height=3.21, dpi=600, units="cm") 




gofile = "/Users/gturco/Documents/code/Brady/sorghum/co-express/slope/association"
y = read.table(gofile, sep="\t", header=FALSE)
z = merge(x, y,  by.x = c("gene"), by.y = c("V1"), all = FALSE)

##9.3 5.56
ylim(15, 20).
l = rbind(y,x,z)
levels(l$Anno) <- c("Cell Wall","Reference", "Hemicellulose","None","TF","Lignin","PCD","Cellulose")
ggplot(data=l, aes(x = tissue, y = value, group=gene, color=Anno)) + geom_line() + theme_classic() + theme(text=element_text(size=7), panel.margin = unit(0, "cm"), plot.margin = unit(c(0.06, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm"),legend.key.height  = unit(.2, "cm")) + xlab(NULL) + ylab("RPKM Mean") + scale_colour_manual(values=c("#8d2b53","black","#d86134", "#626362","#27295d","#00492d","#8d2b53","#872328"))
ggsave("/Users/gturco/Documents/Data/Sorg/legend_coexprss_go.pdf", width=5.60, height=3.21, dpi=600, units="cm") 


###2.205 by 1.268
#### Supplemental Data for dist. of low expressed rpkm values
qplot(rpkm,data=v,geom="histogram") +  theme_classic() + theme()  + theme(text=element_text(size=8), legend.key.height  = unit(.15, "cm"), legend.key.width  = unit(.2, "cm"))
ggsave("/Users/gturco/Documents/Data/Sorg/vas_lowrpkm_dist.pdf", width=4.20, height=2.4, dpi=600, units="cm")

qplot(rpkm,data=n,geom="histogram") +  theme_classic() + theme()  + theme(text=element_text(size=8), legend.key.height  = unit(.15, "cm"), legend.key.width  = unit(.2, "cm"))
ggsave("/Users/gturco/Documents/Data/Sorg/non_lowrpkm_dist.pdf", width=4.20, height=2.4, dpi=600, units="cm")

qplot(rpkm,data=r,geom="histogram") +  theme_classic() + theme()  + theme(text=element_text(size=8), legend.key.height  = unit(.15, "cm"), legend.key.width  = unit(.2, "cm"))
ggsave("/Users/gturco/Documents/Data/Sorg/root_lowrpkm_dist.pdf", width=4.20, height=2.4, dpi=600, units="cm")

qplot(rpkm,data=s,geom="histogram") +  theme_classic() + theme()  + theme(text=element_text(size=8), legend.key.height  = unit(.15, "cm"), legend.key.width  = unit(.2, "cm"))
ggsave("/Users/gturco/Documents/Data/Sorg/shoot_lowrpkm_dist.pdf", width=4.20, height=2.4, dpi=600, units="cm") + scale_colour_manual(values=c("black","#3c4a96","grey","#71b436"))


#### heatmap of meth and genes
file_name = "/Users/gturco//Documents/code/Brady/sorg_figure8/cluster_meth_CHG_1.txt"
x = read.table(file_name ,sep='\t',header=TRUE)
rownames(x) <- x[,1]
colnames(x) <- c("pos","sig")
pheatmap(x[,2:5], color=colorRampPalette(brewer.pal(9,"Blues"))(100))
pheatmap(x[,2:5], color=colorRampPalette(brewer.pal(5,"RdBu"))(100))


file_name = "/Users/gturco/Documents/code/Brady/bs_seq_analysis/genic_analysis/CG_all_2_anova.txt"
x = read.table(file_name ,sep='\t',header=TRUE)
colnames(x) <- c("pos","sig")
x$pval <- cut(x$sig, breaks = c(-Inf, 0.001, 0.005, 0.05, 1))
ggplot(x, aes(pos, "", fill = pval)) + geom_tile() + scale_fill_brewer(palette = "Greens",direction=-1)  +
  theme_classic() + theme(text=element_text(size=20))    +  
  theme(text=element_text(size=7), panel.margin = unit(0, "cm"), 
        plot.margin = unit(c(0, 0, 0.01, -0.1), "cm"),axis.ticks = element_line(size = 0.1), 
        axis.line = element_line(size = 0.1), axis.ticks.length = unit(0.05, "cm"), axis.ticks.margin =unit(0, "cm")) + 
  xlab(NULL) + ylab(NULL) 