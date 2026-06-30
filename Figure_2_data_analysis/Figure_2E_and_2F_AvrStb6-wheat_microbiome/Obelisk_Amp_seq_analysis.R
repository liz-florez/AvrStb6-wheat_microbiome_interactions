#### Metabarcoding amplicon sequence analysis for Florez, et al., 2026 ####
#Obelisk data, 4-8dpi, infected: Mock, IPO323 and IPO323 ΔAvrStb6 

#Packages
library(vegan) 
library(ggplot2)
library(Rmisc)
library(FSA)

#Set PARAMETERS
amplicon=c("16S","ITS2")[1]
method=c("clust")[1]
data=c("obelisk")[1]
sub=c("t4","t8")[1]

#Set DIRECTORIES
result.dir=paste("~/Documents/your_data_folder","your_results_dir",amplicon,sep = "/") ; dir.create(result.dir)
plot.dir=paste(result.dir,"your_plot_dir",sep = "/");dir.create(plot.dir)

#Graphics
colores=c("#56B4E9","#009E73","#D55E00")
tema=theme(axis.text.x = element_text(color="black",size=10, angle=0,family = "Arial"),
           axis.text.y = element_text(color="black",size=10,family = "Arial"),
           axis.title = element_text(color="black",size=10,family = "Arial"),
           legend.text = element_text(color = "black",size=10,family = "Arial"),
           legend.key.size = unit(0.3,"cm"),
           legend.title = element_text(color="black",size=10,family = "Arial"),
           strip.text.x = element_text(size=8, color="black",face="bold",family = "Arial"),
           strip.text.y = element_text(size=8, color="black",face="bold",family = "Arial"),
           legend.position = "right", 
           title = element_text(color="black",size=8,family = "Arial"))
lab=c("Mock", "IPO323", "IPO323 ΔAvrStb6" )

#Load METADATA
met=read.delim(file=paste(result.dir,paste("metadata.table", data,"txt",sep = "."),sep = "/"))
rownames(met)=met$sample.id

met$treatment=factor(met$treatment,levels = c("Mock","IPO323","koAVRStb6"))
met$host.id=factor(met$host.id,levels = c("OB"))
met$time=factor(met$time, levels = c("t4","t8"))

#Load DATA
if (method == "clust" | method == "denoise"){
  if (amplicon=="16S"){
    #16s
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,"txt",sep = "."),sep = "/"))
  } else if (amplicon=="ITS2"){
    #ITS
    otu=read.delim(file=paste(result.dir,paste("otu.table",data,"txt",sep = "."),sep = "/"))
    colnames(otu)=gsub("X","",colnames(otu))
    tax=read.delim(paste(result.dir,paste("taxa.table",data,,"txt",sep = "."),sep = "/"))
  }
}

#Sanity CHECK
sum(rowSums(otu)==0)
sum(colSums(otu)==0)

#Make SUBSET
met.2=met[met$host.id %in% c("OB") & met$treatment %in% c("Mock","IPO323","koAVRStb6") & met$time%in%sub,]
otu.2=otu[,colnames(otu) %in% met.2$sample.id]
met.2=met.2[rownames(met.2) %in% colnames(otu.2),]

otu.2=otu.2[,colSums(otu.2)>100]
otu.2=otu.2[rowSums(otu.2)>0,]
met.2=met.2[colnames(otu.2),]

tax=tax[rownames(otu.2),]
  
# NORMALIZATION for diversity analysis
set.seed(23171341)
otu.s=t(rrarefy(t(otu.2), sample = min(colSums(otu.2))))
print(min(colSums(otu.2)))
otu.r=t(t(otu.s)/colSums(otu.s)*100)

#Reported SUMMARY  
sum(colSums(otu.2))
dim(otu.2)
summary(colSums(otu.2))


#Calculate ALPHA-diversity
ric=colSums(otu.s>0)
shan=diversity(t(otu.s),index = "shannon")

#Plot ALPHA-diversity
alpha=data.frame(richness=ric, shannon=shan, total=colSums(otu.2) )
sum(rownames(alpha)==rownames(met.2))==dim(alpha)[1]
alpha=cbind(alpha, met.2[,c("plate.number","plate.id","well","sample.id","host",
                            "sample.name", "host.id","treatment","time","individual")])

a=ggplot(alpha, aes(x=treatment, y=shannon, fill=treatment))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(alpha=0.50,shape=21,size=1.5,color="black",width = 0.2, aes(fill=treatment))+
  scale_fill_manual(values=colores, labels=lab)+ylab("Shannon index")+xlab("Treatment")+
  theme_minimal()+tema+ggtitle(amplicon)+
  scale_x_discrete(labels=lab)
a

ggsave(paste(sub,method,amplicon,"shannon.treatment2.png",sep = "."),
       path = plot.dir,
       width = 5, height = 2.5, dpi = 300)

#Test ALPHA-diversity
kruskal.test(shannon~treatment, data = alpha)
res=dunnTest(shannon~treatment, data = alpha, method = "bh")$res
res
write.table(res, paste(plot.dir,paste(sub,amplicon,method,"kruskall-Dunn_shannon",sep="."),sep="/"),sep = "\t")

#Calculate NMDS
set.seed(23171341)
scaling=vegdist(t(otu.s), method = "bray", binary = F)
scaling2=metaMDS(scaling) ; scaling2$stress
scaling3=data.frame(scaling2$points) 
scaling3=cbind(scaling3,met.2[rownames(scaling3),])
scaling3=cbind(scaling3,alpha[rownames(scaling3),c("richness","shannon","total")])

#Plot NMDS
a=ggplot(data=scaling3[,], aes(x=MDS1, y=MDS2, shape=host.id,color=treatment, fill=treatment))+
  geom_point(alpha=1,size=3)+theme_minimal()+tema+
  scale_fill_manual(values=c(colores),labels=lab)+ylab("NMDS2")+xlab("NMDS1")+
  scale_color_manual(values=c(colores),labels=lab)+
  scale_shape_manual(values=c(21:25,21))+
  ggtitle(amplicon, "NMDS-Bray")
a

ggsave(paste(sub,method,amplicon,"NMDS2.png",sep = "."),
       path = plot.dir,
       width = 4, height = 2.5, dpi = 300)

#Test BETA-diversity (PERMANOVA)
set.seed(23171341)
scaling=vegdist(t(otu.s), method = "bray", binary = F) 
per=adonis2(scaling~treatment, data = met.2, permutations = 1000)
per
write.table(per, paste(plot.dir,paste(sub,amplicon,method,"permanova.host.treatment.txt",sep="."),sep="/"),sep = "\t")


#Barplot - Relative ABUNDANCE

#Packages
detach("package:Rmisc", unload = TRUE)
detach("package:plyr", unload = TRUE)
library(dplyr)
library(tidyr)

#Prepare data
t="genus"
otu.s2=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s2, key = "sample", value = "abs", colnames(otu.2))
otu.t[is.na(otu.t[,t]),t]="Unclassified"
for (i in met.2$sample.id){
  otu.t[otu.t$sample == i, c( "host.id","time","treatment","plate.id")] = met.2[i,c("host.id","time","treatment","plate.id" )]
}

#Calculate relative abundance by genus
otu.t0 = otu.t %>% group_by(host.id,treatment,time) %>% summarise(absab = sum(abs))
otu.t1 = otu.t %>% group_by(otu.t[,t],host.id,treatment,time) %>% summarise(absab = sum(abs)) 
colnames(otu.t1)[1]=c("rank")

sum(otu.t1$host.id==otu.t0$host.id, na.rm = T)==dim(otu.t1)[1] #Sanity check
sum(paste(otu.t1$host.id,otu.t1$rank)==paste(otu.t0$host.id,otu.t1$rank))==dim(otu.t1)[1] #Sanity check

otu.t1$relab=otu.t1$absab/otu.t0$absab*100
otu.t1$rank=gsub("g:","",otu.t1$rank)

#Open colors and isolate genus 
col=read.delim(paste(result.dir,"genera_colors",sep = "/"))

#Reorder taxa if needed
orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
high = orden[orden$relab >= 1, "rank"] #filter high bundance otus
p=sort(unique(c(high$rank,col$rank)))
p=p[!p%in%c("Gp6","Gp16")]

otu.t1[!otu.t1$rank %in% p, "rank" ] = "Other"
niveles=c(sort(p[!(p%in% c("Unclassified", "Other"))]),"Unclassified","Other")
otu.t1$rank=factor(otu.t1$rank, levels = niveles)

a=ggplot(data=otu.t1, aes(y=relab, x=treatment, fill=rank))+
  geom_bar(stat="identity",linewidth = 0.95,size=0.5)+
  scale_fill_manual(values = col[col$rank %in% unique(otu.t1$rank),"color"] )+
  ylab("Relative abundance (%)")+ xlab("Treatment")+labs(fill="Genus")+
  #facet_wrap(~host.id*time, scales="free_x", ncol = 4)+
  theme_minimal()+tema+
  scale_x_discrete(labels=lab)

a

ggsave(paste(sub,method,amplicon,t,"barplot.all.custom_new_theme_newcolors2.png",sep = "."),
       path = plot.dir,
       width = 5.5, height = 2.5, dpi = 300)

# DIFERENTIAL abundance analysis

library(metagenomeSeq)
library(phyloseq)

#Settings 
colores3=col$color ; names(colores3)=col$rank
time=sub
host=c("OB")
treat=c("Mock","IPO323","koAVRStb6")
cut=0.5 #frequency (prevalence) cut off 

#Create phyloseq object and set levels
g=phyloseq(otu_table(as.matrix(otu.2), taxa_are_rows = T), sample_data(alpha), tax_table(as.matrix(tax)))
sample_data(g)$treatment = factor(sample_data(g)$treatment, levels = treat)
sample_data(g)$time = factor(sample_data(g)$time,levels = time)
sample_data(g)$host.id = factor(sample_data(g)$host.id,levels = host)
ranks=c("domain","phylum","class","order","family","genus")

#Set results directory
dir.1=paste(plot.dir,"enrichments",sep = "/") #name for the main folder
dir.create(dir.1)

dir.2=paste(dir.1,"between_treatments_final",sep = "/") #Dont change this name
dir.create(dir.2)

#Test Mock vs Mutant

#Subset phyloseq objects
g1=subset_samples(g, treatment %in% c("Mock","koAVRStb6"))
g1=prune_taxa((taxa_sums(g1)>2), g1)
counts=as( otu_table( g1 ), "matrix" )
taxon=AnnotatedDataFrame(as.data.frame(tax_table(g1)))
experim=AnnotatedDataFrame(as.data.frame(sample_data(g1)))

#Normalization 
mp=newMRexperiment(counts = counts, phenoData = experim, featureData = taxon)
mp.n=cumNorm( mp, p=cumNormStat( mp ) )

#DAA gausian model
t = pData(mp.n)$treatment
normFactor = normFactors(mp.n)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~t + normFactor)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = mp.n, mod = mod, useCSSoffset = FALSE, control = settings)

#Summarise result
res=MRfulltable(fit, number = dim(MRcounts(mp.n))[1])
res$samples=res$`+samples in group 0` + res$`+samples in group 1`
res2=res[res$samples>=dim(sample_data(g1))[1]*cut & res$adjPvalues<0.05 & abs(res$tkoAVRStb6)>1,]
res2=cbind(res2,tax[rownames(res2),])

write.table(res, paste(dir.2,paste(sub,"fulltable","Mock_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)

write.table(res2, paste(dir.2,paste(sub, "cutofftable","Mock_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)

#Plot DAA
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"family"] #Quick fix for unclassified genera
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"order"]
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"class"]
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"phylum"]
res2$genus=gsub("_genera_incertae_sedis","",res2$genus)

res2$genus=gsub("[a-z]:","",res2$genus)
res2$rank=res2$genus
res2$rank=paste(res2$genus,res2$otu.id,sep = ":")
res2$rank=factor(res2$rank, levels = res2[order(res2$tkoAVRStb6,decreasing = F),"rank"])

a=ggplot(data=res2, aes(x=(tkoAVRStb6)*-1, y=rank, color=genus, size=samples/dim(sample_data(g1))[1]*100))+
  geom_point()+scale_color_manual(values=colores3 )+
  geom_vline(xintercept = 0)+scale_size_continuous(limits=c(50,100))+
  scale_x_continuous(limits = c(-7.5,7.5))+scale_y_discrete(position = "left")+
  theme_minimal() + tema +
  labs(x = expression(log[2]*" fold change "),
       y = "Differential OTU",
       color = "Genus",
       size = "Prevalence (%)", 
       title = "IPO323 ΔAvrStb6 vs Mock"
       ) +
  scale_x_continuous(breaks = c(-6,-4,-2,0,2,4,6), limits = c(-7.5,7.5))

a

ggsave(paste(plot.dir,
             paste(sub, "cutofftable","Mock_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.png",sep="_"),
             sep = "/") , width = 5.5, height = 2.6, dpi = 300)


#Test Mock vs Wt

#Subset phyloseq objects
g1=subset_samples(g, treatment %in% c("Mock","IPO323"))
g1=prune_taxa((taxa_sums(g1)>3), g1)
counts=as( otu_table( g1 ), "matrix" )
taxon=AnnotatedDataFrame(as.data.frame(tax_table(g1)))
experim=AnnotatedDataFrame(as.data.frame(sample_data(g1)))

#Normalization 
mp=newMRexperiment(counts = counts, phenoData = experim, featureData = taxon)
mp.n=cumNorm( mp, p=cumNormStat( mp ) )
#mat = MRcounts(mp.n, norm = TRUE, log = TRUE)

#DAA gausian model
t = pData(mp.n)$treatment
normFactor = normFactors(mp.n)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~t + normFactor)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = mp.n, mod = mod, useCSSoffset = FALSE, control = settings)

#Summarise results
res=MRfulltable(fit, number = dim(MRcounts(mp.n))[1])
res$samples=res$`+samples in group 0` + res$`+samples in group 1`
res2=res[res$samples>=dim(sample_data(g1))[1]*cut & res$adjPvalues<0.05 & abs(res$tIPO323)>1,]
res2=cbind(res2,tax[rownames(res2),])

write.table(res, paste(dir.2,paste(sub,"fulltable","Mock_vs_IP0323",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)

write.table(res2, paste(dir.2,paste(sub,"cutofftable","Mock_vs_IP0323",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)

#PLot DAA
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"family"] #Quick fix for unclassified genera
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"family"]
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"order"]
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"class"]
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"phylum"]

res2$genus=gsub("[a-z]:","",res2$genus)
res2$rank=res2$genus
res2$rank=paste(res2$genus,res2$otu.id,sep = ":")
res2$rank=factor(res2$rank, levels = res2[order(res2$tIPO323,decreasing = F),"rank"])

a=ggplot(data=res2, aes(x=(tIPO323)*-1, y=rank, color=genus, size=samples/dim(sample_data(g1))[1]*100))+
  geom_point()+scale_color_manual(values=colores3 )+
  geom_vline(xintercept = 0)+scale_size_continuous(limits=c(50,100))+
  scale_x_continuous(limits = c(-7.5,7.5))+scale_y_discrete(position = "left")+
  theme_minimal() + tema +
  labs(x = expression(log[2]*" fold change "),
       y = "Differential OTU",
       color = "Genus",
       size = "Prevalence (%)", 
       title = "IPO323 vs Mock") +
  scale_x_continuous(breaks = c(-6,-4,-2,0,2,4,6), limits = c(-7.5,7.5))

a

ggsave(paste(plot.dir,
             paste(sub, "cutofftable","Mock_vs_IPO323",cut,"OTU.ID","CSS_gausian.png",sep="_"),
             sep = "/") , width = 5.5, height = 2.6, dpi = 300)


#Test Wt vs Mutant

#Subset phyloseq objects
g1=subset_samples(g, treatment %in% c("IPO323","koAVRStb6"))
g1=prune_taxa((taxa_sums(g1)>2), g1)

counts=as( otu_table( g1 ), "matrix" )
taxon=AnnotatedDataFrame(as.data.frame(tax_table(g1)))
experim=AnnotatedDataFrame(as.data.frame(sample_data(g1)))

#Normalization
mp=newMRexperiment(counts = counts, phenoData = experim, featureData = taxon)
mp.n=cumNorm( mp, p=cumNormStat( mp ) )

#DAA gausian model
t = pData(mp.n)$treatment
normFactor = normFactors(mp.n)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~t + normFactor)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = mp.n, mod = mod, useCSSoffset = FALSE, control = settings)

#Summarise results
res=MRfulltable(fit, number = dim(MRcounts(mp.n))[1])
res$samples=res$`+samples in group 0` + res$`+samples in group 1`
res2=res[res$samples>=dim(sample_data(g1))[1]*cut & res$adjPvalues<0.05 & abs(res$tkoAVRStb6)>1,]
res2=cbind(res2,tax[rownames(res2),])

#Plot DAA
res2[is.na(res2$genus),"genus"]=res2[is.na(res2$genus),"family"] #Quick fix for unclassified genera
res2$genus=gsub("[a-z]:","",res2$genus)
res2$rank=res2$genus
res2$rank=paste(res2$genus,res2$otu.id,sep = ":")
res2$rank=factor(res2$rank, levels = res2[order(res2$tkoAVRStb6,decreasing = F),"rank"])

a=ggplot(data=res2, aes(x=(tkoAVRStb6)*-1, y=rank, color=genus, size=samples/dim(sample_data(g1))[1]*100))+
  geom_point()+scale_color_manual(values=colores3 )+
  geom_vline(xintercept = 0)+scale_size_continuous(limits=c(50,100))+
  scale_x_continuous(limits = c(-7.5,7.5))+scale_y_discrete(position = "left")+
  theme_minimal() + tema +
  labs(x = expression(log[2]*" fold change "),
       y = "Differential OTU",
       color = "Genus",
       size = "Prevalence (%)", 
       title = "IPO323 ΔAvrStb6 vs IPO323")+
  scale_x_continuous(breaks = c(-6,-4,-2,0,2,4,6), limits = c(-7.5,7.5))

a

ggsave(paste(plot.dir,
             paste(sub, "cutofftable","IPO323_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.png",sep="_"),
             sep = "/") , width = 5.5, height = 2.6, dpi = 300)


write.table(res, paste(dir.2,paste(sub,"fulltable","IP0323_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)

write.table(res2, paste(dir.2,paste(sub,"cutofftable","IPO323_vs_koAVRStb6",cut,"OTU.ID","CSS_gausian.txt",sep = "_"),sep = "/"),
            quote = F, sep = "\t", row.names = F, col.names = T)


###############################

