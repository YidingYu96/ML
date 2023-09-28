#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

#引用包
library(ggpubr)          
library(limma)
library(ggsci)

#输入文件  
expfile="GSE1.txt"          #表达矩阵
hub="gene.txt"                     #核心基因
sample="sample.txt"               #样本属性
Ccol="blue"                       #正常控制组颜色
Pcol="red"                        #疾病实验组颜色
afmethod="t.test"        #t检验："t.test"（小样本量）,秩和检验:"wilcox.test"（大样本量）

rt=read.table(expfile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=t(rt)
rt2=read.table(sample,sep="\t",header=F,check.names=F,row.names = 1)
ssample=as.character(intersect(rownames(rt2),rownames(rt)))
rt=rt[ssample,]
rt2=rt2[ssample,,drop=F]
rt1=cbind(rt2,rt)
colnames(rt1)[1]="Type"
rt1$ID=rownames(rt1)
genes=read.table(hub,sep="\t",header=F,check.names=F)[,1]
rt1=rt1[,c("ID","Type",genes)]
af=colnames(rt1)
df=c(3:ncol(rt1))

for (i in df) {

outFile=paste0(af[i],"_","boxplotdiff.pdf")             
rt=as.data.frame(rt1[,c(1,2,i)])
x=colnames(rt)[2]
y=colnames(rt)[3]
colnames(rt)=c("id","Type","Expression")


group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
rt$Expression=as.numeric(rt$Expression)

comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


boxplot=ggboxplot(rt, x="Type", y="Expression", fill="Type",
		          xlab=x,     
		          ylab=y,
		          legend.title=x,
		          palette = c(Ccol,Pcol),   
		          #add = "jitter"
		          )+ 
	    stat_compare_means(comparisons = my_comparisons,method=afmethod)


pdf(file=outFile,width=4,height=4)
print(boxplot)
dev.off()

}

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
