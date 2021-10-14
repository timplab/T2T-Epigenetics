library(tidyverse)
library(mclust)

gene.meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/chr7/CGI_GeneBody_SingleRead_methylation1kbGB.tsv")

df.cgi <- gene.meth %>%
  dplyr::select(c(transcript_id,qname,cgi)) %>%
  unite(newname,transcript_id,qname, sep=":") 

class <- df.cgi$newname
X <- df.cgi[,-1]
model1 <- Mclust(X,G=2)
summary(model1)
#plot(model1, what = c("classification"))
df.cgi$cgi_cluster <- model1$classification



df.gb <- gene.meth %>%
  select(c(transcript_id,qname,gb))  %>%
  unite(newname1,transcript_id,qname, sep=":") 

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

write.table(clustered_reads, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/revision_analysis/promoter_clustering/ClusteredReads_chr7_proteinCoding.tsv",
            quote=F,
            row.names=F,
            sep="\t")

ggplot(clustered_reads,aes(x=cgi,y=gb,color=as.factor(quantile)))+geom_point(alpha=.05)

stat.sum <- clustered_reads %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(sum=n())

stat <- clustered_reads %>%
  group_by(name) %>%
  mutate(total=n()) %>%
  ungroup() %>%
  group_by(name,group,quantile,total) %>%
  summarise(sum=n()) %>%
  mutate(frac=sum/total) %>%
  distinct()


ggplot(stat, aes(x=as.factor(quantile),y=frac,color=as.factor(quantile)))+geom_violin(width=1)+geom_boxplot(width=.1,outlier.shape = NA)+facet_wrap(~group)+theme_bw() 

