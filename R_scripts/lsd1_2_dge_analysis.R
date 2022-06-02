################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: lsd1_2_dge_analysis.R
# 
# 
# 
# input: gene quantification file generated from  featureCounts_gene.sh
# output:  DGE output file
# required software: R version 4.0.4
################################################################################
library(edgeR)

#functions here
readit = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("Aligned.sortedByCoord.out.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}


dgplus = function(df,group1,group2,grouping){
  grp1 = group1
  grp2 = group2
  #subset original dataframe by merging string groups
  temp1 = subset(df,select=grp1)
  temp2 = subset(df,select=grp2)
  #new subset dataframe below
  df = transform(merge(temp1,temp2,by="row.names"),row.names=Row.names,Row.names=NULL)
  
  y <- df
  y_full <- DGEList(counts=y,group = grouping)
  keep <- rowSums(cpm(y_full) > 0) >= 3
  lm1.y <- y_full[keep, keep.lib.sizes =FALSE]
  lm1.y <- calcNormFactors(lm1.y, method="TMM")
  lm1.tmm = as.data.frame(lm1.y$counts)
  lm1.design <- model.matrix( ~ as.numeric(group==1), data = lm1.y$samples)
  lm1.v <- voom(lm1.y, lm1.design, plot=FALSE)
  lm1.fit <- lmFit(lm1.v, lm1.design)
  lm1.fit.c <- eBayes(lm1.fit)
  lm1.top <- topTable(lm1.fit.c, sort="P",number='all')
  lm1.top.merge = transform(merge(lm1.top,df,by="row.names"),row.names=Row.names,Row.names=NULL)
  lm1.top.merge = lm1.top.merge[order(lm1.top.merge$adj.P.Val),]
  return(lm1.top.merge)
}


pombe_geneOfrac = readit("Spombe_gene_withO_fraction.txt")
#sample labels
controls = c("C1","C2","C3")
lsd1c1 = c("T7","T8","T9")
lsd1c2 = c("T10","T11","T12")
#grouping for design matrix
lsd1c1_grp = ifelse(c("C1","C2","C3","T7","T8","T9") == lsd1c1,1,0)
lsd1c2_grp = ifelse(c("C1","C2","C3","T10","T11","T12") == lsd1c2,1,0)
lsd1.geneOfraction.plus = dgplus(pombe_geneOfrac,controls,lsd1c1,lsd1c1_grp)
lsd2.geneOfraction.plus = dgplus(pombe_geneOfrac,controls,lsd1c2,lsd1c2_grp)
#output DGE lists for each mutant type
write.table(lsd1.geneOfraction.plus,file="lsd1_dge.txt",quote=F,sep="\t",row.names=T,col.names=T)
write.table(lsd2.geneOfraction.plus,file="lsd2_dge.txt",quote=F,sep="\t",row.names=T,col.names=T)

