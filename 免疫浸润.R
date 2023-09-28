# install packages 这三个安装不成功的话，就安后面的bseqsc包也行
#install.packages('e1071')
#install.pacakges('parallel')
#install.packages('preprocessCore')
library(e1071)
library(preprocessCore)
library(parallel)

#install.packages('devtools')
library(devtools)
#devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)#这个包携带大量CIBERSORT的依赖，前三个安装不好可以安装他

################安装CIBERSORT包##########################################################
#if(!require(CIBERSORT))devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
# 包全部安装完成
#BiocManager::install("ComplexHeatmap")
# 画热图的包
#install.packages("pheatmap")
#install.packages("ComplexHeatmap")
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
# 同时准备好LM22的TXT文件，注意自己以后的文件要和这个TXT的格式一样
# 加载CIBERSORT包成功后，系统内部会自带data(LM22)
data(LM22) 
#data(GSE1)#TCGA的演示数据，正式情况下就用自己的数据
inputFile="GSE1.txt"
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rownames(rt) = rt[,1]
rt = rt[,-1]
# 正式开始探索
# 看5*5的数据
LM22[1:5,1:5]
rt[1:5,1:5]

# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
results <- cibersort(sig_matrix = LM22, mixture_file = rt,perm = 1000, QN = F)
write.table(results,"results.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# 理解一下results的结果
# 你可以理解为返回了一个列名为细胞类型、行名为样本名的细胞浸润程度（占比）的矩阵
# 此外result中还会多出三列：
# P-value: 用来展示去卷积的结果在所有细胞类群中是否具有差异
# Correlation:参考矩阵与输入矩阵的特征基因相关性
# RMSE: Root mean squared error，参考矩阵与输入矩阵的特征基因标准差

# heatmap
# 按行（样本内部）标准化可以看出在各类样本内部，M2浸润程度（占比）最高
rowscale <- results[,1:ncol(LM22)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=F,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度

# 各类样本之间也具有自己占比高的特异性免疫细胞
columnscale <- results[,1:ncol(LM22)]
columnscale <- columnscale[,apply(columnscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(columnscale,
         scale = 'column',
         cluster_col=F,
         cluster_row=T,
         angle_col = "315")

# 堆积比例图
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175'
)
cellnum <- results[,1:ncol(LM22)]
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],rownames(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}

colnames(data4plot)<-c('proportion','celltype','sample')
data4plot$proportion <- as.numeric(data4plot$proportion)
ggplot(data4plot,aes(sample,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my36colors)+#自定义fill的颜色
  ggtitle("cell portation")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x=element_text(size=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+#把x坐标轴横过来
  guides(fill=guide_legend(title=NULL))


#t检验

#第一列是iD,group在第二列

library(tidyverse)
library(ggpubr)
library(dplyr)
library(rstatix)
library(reshape2)
exp <- read.table("ttestresult.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

mydata.long <- exp %>%
  pivot_longer(-group , names_to = "variables", values_to = "value")

stat.test <- mydata.long%>%
  group_by(variables)%>%
  t_test(value ~ group)%>%
  adjust_pvalue(method ="fdr")%>%
  add_significance()
stat.test

write.table(stat.test, file = "tresult.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

