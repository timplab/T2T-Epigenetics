library(tidyverse)
library(mclust)

gene.meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/RNAseq/CGI_GeneBody_SingleRead_methylationRNAseqchr7_HG002.tsv")

df.cgi <- gene.meth %>%
  select(c(name,qname,cgi)) %>%
  unite(newname,name,qname, sep=":") 

class <- df.cgi$newname
X <- df.cgi[,-1]
model1 <- Mclust(X,G=2)
summary(model1)
#plot(model1, what = c("classification"))
df.cgi$cgi_cluster <- model1$classification



df.gb <- gene.meth %>%
  select(c(name,qname,gb))  %>%
  unite(newname1,name,qname, sep=":") 

class <- df.gb$newname1
X <- df.gb[,-1]
model1 <- Mclust(X,G=2)
summary(model1)
plot(model1, what = c("classification"))
df.gb$gb_cluster <- model1$classification


ggplot(df.gb,aes(x=gb,fill=as.factor(gb_cluster)))+geom_histogram()
ggplot(df.cgi,aes(x=cgi,fill=as.factor(cgi_cluster)))+geom_histogram()

clustered_reads <- cbind(df.gb,df.cgi) %>%
  filter(newname == newname1) %>%
  separate(newname, into=c("name","qname"),":") %>%
  mutate(cgi_cluster=ifelse(cgi_cluster==1, "low_cgi","high_cgi")) %>%
  mutate(gb_cluster=ifelse(gb_cluster==1, "low_gb","high_gb")) %>%
  unite(group,cgi_cluster,gb_cluster, sep=":") 

clustered_reads$quantile = gene.meth$quantile

write.table(clustered_reads, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/chr7/ClusteredReads_chr7_HG002_proteinCoding_RNAseq.tsv",
            quote=F,
            row.names=F,
            sep="\t")

ggplot(clustered_reads,aes(x=cgi,y=gb,color=as.factor(quantile)))+geom_point(alpha=.05)

stat.sum <- clustered_reads %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(sum=n()) 

id <- clustered_reads %>%
  select(c(name,quantile)) %>%
  distinct()

stat <- clustered_reads %>%
  group_by(name) %>%
  mutate(total=n()) %>%
  group_by(name,group,quantile,total) %>%
  summarise(sum=n()) %>%
  ungroup() %>%
  select(c(name,group,total,sum)) %>%
  complete(name,group,fill=list(sum=0)) %>%
  mutate(sum=ifelse(sum == 1, 0, sum)) %>%
  mutate(frac=sum/total) %>%
  filter(total > 10) %>%
  distinct()

stat[is.na(stat)] <- 0
stat <- merge(stat,id, by="name")

ggplot(stat, aes(x=as.factor(quantile),y=frac,color=as.factor(quantile)))+geom_violin(width=1)+geom_boxplot(width=.1,outlier.shape = NA)+facet_wrap(~group)+theme_bw()
