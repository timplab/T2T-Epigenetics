asm=c("chm13", "grch38", "chm13", "grch38")
sample=c("chm13", "chm13", "hg002", "hg002")
num=c(32186696,
      28979828,
      29075692,
      28259048)
seq=c("ont", "ont", "ill", "ill")

df = data.frame(asm=asm, sample=sample, seq=seq,num=num)


p1 <- ggplot(data=df, aes(y=num/1e6, x=asm,group=seq, color=asm,shape=seq))+geom_point()+geom_path(position = "identity")+ theme(text = element_text(size=20))+theme_classic()+labs( y = "Number of MACs Peaks")+theme(legend.position = "bottom")


ggsave(filename = paste0(figs, "/Number_CpGs.pdf"), plot = p1,
       width = 5,
       height = 5)