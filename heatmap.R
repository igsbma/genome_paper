source("/Users/bma/Dropbox (IGS)/Projects/_technical/R_Library/heatmapLib2.R")
source("/Users/bma/Dropbox (IGS)/Projects/_technical/R_Library/vegdist.R")

################# set 3 ###############
mytotal=read.csv("set3.csv",header=T,row.names=1,check.names=F)
dim(mytotal)
head(mytotal)
mytotal[1:5,1:5]

mydata=mytotal
mydata=mytotal[,-2]
mydata[1:5,1:5]

groupf=factor(CO$group)
levels(groupf)
groupTbl=c()
groupTbl[levels(groupf)[1]]="navy"
groupTbl[levels(groupf)[2]]="green"
groupTbl[levels(groupf)[3]]="red"
groupTbl[levels(groupf)[4]]="yellow"

cs=colSums(mydata)
o=order(cs,decreasing=T)
ct=mydata[,o]
ct[1:5,1:5] 
pt=t(apply(ct,1,function(x) 100*x/sum(x)))
pt[1:5,1:5]
nCols=25
dim(pt)
hc.ward=hclust(dist(pt),method="ward.D")

nClrs=3
memb.ward=cutree(hc.ward,k=nClrs)
table(memb.ward)

lx=0.86
sideBars.ward=cbind(colorTbl[memb.ward],groupTbl[groupf])
colnames(sideBars.ward)=c("Cluster","group")
rownames(sideBars.ward)=memb.ward

write.csv(file="cluster_ward.csv", memb.ward)

nCols=23
dim(pt)
postscript("./heatmap_ward.eps",hor=T)

heatmap2(as.matrix(pt)[,1:nCols],         # show only nCols first columns of the table tbl
         col=rainbow(50,start=1/6,end=0),
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         Rowv=as.dendrogram(hc.ward),  
         RowSideColors=sideBars.ward,    # RowSideColors=sideBars,          # this sets three color columns: clustering, pH and Nugent score
         RowSideTitleCex=0.5,             # scaling factor for the titles of the color bars
         RowSideTitleLine=0.5,            # the distance of the titles for the color bars
         margins=c(18,4),                 # =c(bottom margin, right margin)
         labRow=NA,                       # suppress row labels
         cexCol=1.6,                      # magnification factor for column labels
         xlas=2,                          # column labels are to be perpendicular to the x-axis
         main="")                         # main titlek#legend(lx,0.97,legend=c("Not Transmitted","Transmitted"),bg="white",fill=transmissionTbl,title="transmission",cex=0.7)
legend(lx,0.97,legend=c("ATCC", "DESULPHO", "MD","PBS"),bg="white",fill=groupTbl,title="Group",cex=0.9)
dev.off()



