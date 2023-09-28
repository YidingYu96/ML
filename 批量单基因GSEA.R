#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
install.packages("devtools")
devtools::install_github("YuLab-SMU/DOSE")
devtools::install_github("YuLab-SMU/HDO.db")   
devtools::install_github('YuLab-SMU/clusterProfiler')


library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)

#输入文件
expFile="GSE1.txt"          #表达矩阵，去除正常控制组
hub="gene.txt"    #核心基因文件名称

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneaf=read.table(hub,sep="\t",header=F,check.names=F)[,1]


for (genei in geneaf) {
  

#按基因表达分组
group <- ifelse(data[c(genei),]> median(data[c(genei),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

#差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
#保存单基因分组的所有基因差异结果
#DIFFOUT=rbind(id=colnames(Diff),Diff)
#write.table(DIFFOUT,file=paste0("1.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)


#展示差异最大的前30个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]
#分组标签
Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(High="red",Low="blue"  ))



logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.",paste0(genei,"_all_GSEA.xls")),sep="\t",quote=F,col.names=T)

#排序后分别取GSEA结果的前5个和后5个
num=5
pdf(paste0("2.",paste0(genei,"_down_GSEA.pdf")),width = 8,height = 8)
af=gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
print(af)
dev.off()
pdf(paste0("2.",paste0(genei,"_up_GSEA.pdf")),width = 8,height = 8)
af=gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
print(af)
dev.off()
#排序后取前5个和后5个一起展示
num=5
pdf(paste0("2.",paste0(genei,"_all_GSEA.pdf")),width = 10,height = 10)
af=gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
print(af)
dev.off()


}

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af