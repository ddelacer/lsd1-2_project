################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: lsd1_2_visuals.R
# 
# 
# 
# input: DGE output file from lsd1_2_dge_analysis.R
# output:  Volcano plots
# required software: R version 4.0.4
################################################################################
library(ggplot2)
library(VennDiagram)
library(extrafontdb)
library(extrafont)
#functions
vplot = function(df, ttext,stext){
  txt1= ttext
  txt2 = stext
  # color="#0743FA" #higlighting higlighting where l1 is >2 or <-2 (blue)
  # color="#FF2600" #higlighting where l2 is >2 or <-2 (red)
  #ggplot(df,aes(logFC,-log10(P.Value),color="#0743FA")) + #blue for lsd1
  #ggplot(df,aes(logFC,-log10(P.Value),color="#FF2600")) + #red for lsd2
  ggplot(df,aes(logFC,-log10(P.Value),color=-log10(P.Value)))+
    geom_point(shape=16,size=4, show.legend=F,color="#FF2600") + 
    #geom_point(data=subset(df,-log10(P.Value)<3) , color="gray",size=2) +
    geom_point(data=subset(df,logFC<=2 & logFC>0) , color="gray",size=4) +
    geom_point(data=subset(df,logFC>=-2 & logFC<0) , color="gray",size=4) +
    #geom_vline(xintercept = c(-2,2))+
    scale_x_continuous(name= "lFC",lim=(c(-5,7.8))) +
    scale_y_continuous(name=expression(-log[10]* " P-value"),lim=c(0,10.5)) +
    ggtitle(bquote(italic(.(txt1))),subtitle= paste(  sum(df[,"logFC"]>=2), "lFC >= 2,", sum(df[,"logFC"]<=-2),"lFC <= -2" )) +
    theme(plot.title = element_text(size=30,face="italic",family="Arial"),
          plot.subtitle = element_text(size=15,family="Arial"),
          axis.title.x = element_text(size=30,family="Arial"),
          axis.line.x = element_line(color="black"),
          axis.text.x = element_text(size=15,family="Arial"),
          axis.title.y = element_text(size=30,family="Arial"),
          axis.line.y = element_line(color="black"),
          axis.text.y = element_text(size=15,family="Arial"),
          panel.background = element_rect(fill="white"),
          panel.grid.minor = element_line("white"),
          panel.grid.major = element_line("white")
    )
  
}
#same volcano plot function, just with stripped text
vstrip = function(df, ttext,stext){
  txt1= ttext
  txt2 = stext
  ggplot(df,aes(logFC,-log10(P.Value),color=-log10(P.Value)))+
    geom_point(shape=16,size=4, show.legend=F,color="#FF2600") + 
    geom_point(data=subset(df,logFC<=2 & logFC>0) , color="gray",size=4) +
    geom_point(data=subset(df,logFC>=-2 & logFC<0) , color="gray",size=4) +
    scale_x_continuous(name= "lFC",lim=(c(-5,7.8))) +
    scale_y_continuous(name=expression(-log[10]* " P-value"),lim=c(0,10.5)) +
    ggtitle(bquote(italic(.(txt1))),subtitle= paste(  sum(df[,"logFC"]>=2), "lFC >= 2,", sum(df[,"logFC"]<=-2),"lFC <= -2" )) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black"),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_line(color="black"),
          axis.text.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_line("white"),
          panel.grid.major = element_line("white")
    )
  
}


#functions for pulling out gene lists
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}
Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

vennlist = function(comparelist){
  combs <- 
    unlist(lapply(1:length(comparelist), 
                  function(j) combn(names(comparelist), j, simplify = FALSE)),
           recursive = FALSE)
  #make sure list is correct
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  #get the names of them?
  #str(combs)
  #put the names in elements
  elements <- 
    lapply(combs, function(i) Setdiff(comparelist[i], comparelist[setdiff(names(comparelist), i)]))
  return(elements)
}

#FC gene overlap between mutant experiments
lsd1vslsd2 = venn.diagram(list(
  "Lsd1 lFC" = rownames(lsd1.geneOfraction.plus[lsd1.geneOfraction.plus$logFC>=2 | lsd1.geneOfraction.plus$logFC<=-2,]),
  "Lsd2 lFC" = rownames(lsd2.geneOfraction.plus[lsd2.geneOfraction.plus$logFC>=2 | lsd2.geneOfraction.plus$logFC<=-2,])
),fill=c("#2E86C1","#A569BD"),
alpha=c(.3,.3),
cex = 3,
lwd = .25,
col = "#34495E",
fontfamily = "Arial",
cat.cex = 3,
cat.fontfamily = "Arial",
cat.pos=c(-22,20),
cat.dist=c(.08,.1),
filename=NULL)
dev.off()
tiff(filename="Comparelsd1andlsd1d_8.22.19.tiff",height = 500,width = 900)
grid.draw(lsd1vslsd2)
dev.off()

lsd1v = vplot(lsd1.geneOfraction.plus,"lsd1","")
lsd2v = vplot(lsd2.geneOfraction.plus,"lsd2","")

lsd1vstrip = vstrip(lsd1.geneOfraction.plus,"lsd1","")
lsd2vstrip = vstrip(lsd2.geneOfraction.plus,"lsd2","")

tiff(filename = "lsd1volcano_9.2.19.tiff",height = 400,width = 600)
lsd1v
dev.off()
tiff(filename = "lsd2volcano_9.2.19.tiff",height = 400,width = 600)
lsd2v
dev.off()
tiff(filename = "lsd1volcano_strip_9.2.19.tiff",height = 400,width = 600)
lsd1vstrip
dev.off()
tiff(filename = "lsd2volcano_strip_9.2.19.tiff",height = 400,width = 600)
lsd2vstrip
dev.off()

#combine mutant DGE analyses
lsd12lfcmerge = transform(merge(lsd1.geneOfraction.plus[lsd1.geneOfraction.plus$P.Value<.001,],lsd2.geneOfraction.plus[lsd2.geneOfraction.plus$P.Value<.001,],by="row.names"),row.names=Row.names,Row.names=NULL)

tiff(filename="FC_lsd12_Comparison_strip_9.3.19.tiff",height = 900,width = 900)
ggplot(lsd12lfcmerge,aes(x,y,color=y))+
  geom_point(shape=16,size=9, show.legend=F) + 
  geom_point(data=lsd12lfcmerge[intersect(rownames(lsd12lfcmerge),l1only),] , color="#0743FA",size=9,alpha=.8) + #higlighting higlighting where l1 is >2 or <-2 (blue)
  geom_point(data=lsd12lfcmerge[intersect(rownames(lsd12lfcmerge),l2only) ,] , color="#FF2600",size=9,alpha=.8) + #higlighting where l2 is >2 or <-2 (purple)
  scale_x_continuous(name=expression( "lsd1 log"[2]~"FC"),lim=(c(-5,7.8))) +
  scale_y_continuous(name=expression("lsd2 log"[2]~"FC"),lim=(c(-5,7.8))) +
  ggtitle(expression(paste(italic("lsd1"),"and ",italic("lsd2"),"single mutant comparisons")),subtitle =paste(length(lsd12lfcmerge[,1]),
                                                                                                              "genes",
                                                                                                              "R-squared",
                                                                                                              round(summary(lm(lsd12lfcmerge[,2]~lsd12lfcmerge[,1]))$r.squared,3),
                                                                                                              "P-Value < 2e-16")) +
  geom_smooth(method="lm",formula=y~x,show.legend = F,linetype="solid",size=.5,alpha=.8,color="black")+
  #geom_abline(intercept=0,slope=1,show.legend = F,size=.5,alpha=.5,color="red")+
  scale_color_gradientn(colors=gray(.8))+
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.minor = element_line("white"),
        panel.grid.major = element_line("white")
  )
dev.off()
