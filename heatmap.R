# this is a commend line
source("/Users/bma/Dropbox (IGS)/Projects/_technical/R_Library/heatmapLib2.R")
source("/Users/bma/Dropbox (IGS)/Projects/_technical/R_Library/heatmapLib_w_scale_bar.R")
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
#pt=t(apply(ct,1,function(x) 100*x/sum(x)))
pt=ct #metaPhlAn result
pt[1:5,1:5]
nCols=25
dim(pt)
#hc.comp=hclust(dist(pt),method="comp")
hc.ward=hclust(dist(pt),method="ward.D")

nClrs=3
memb.ward=cutree(hc.ward,k=nClrs)
table(memb.ward)

lx=0.86
sideBars.ward=cbind(colorTbl[memb.ward],siteTbl[sitef],groupTbl[groupf])
colnames(sideBars.ward)=c("Cluster","site","group")
rownames(sideBars.ward)=memb.ward

write.csv(file="cluster_ward.csv", memb.ward)

nCols=15
dim(pt)
postscript("./set3/heatmap_ward_CO2.eps",hor=T)

heatmap2(as.matrix(pt)[,1:nCols],         # show only nCols first columns of the table tbl
         col=rainbow(50,start=1/6,end=0),
         #col=greenred(100), # set heatmap color palette
         #col=colorRampPalette(c("darkblue","black","yellow"))(100),
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         #Rowv = NA,
         Rowv=as.dendrogram(hc.ward),  
         #Rowv=as.dendrogram(hc.ward),    # use our h. clustering
         RowSideColors=sideBars.ward,    # RowSideColors=sideBars,          # this sets three color columns: clustering, pH and Nugent score
         RowSideTitleCex=0.5,             # scaling factor for the titles of the color bars
         RowSideTitleLine=0.5,            # the distance of the titles for the color bars
         margins=c(18,4),                 # =c(bottom margin, right margin)
         #labRow=rownames(mydata),         # add row labels
         labRow=NA,                       # suppress row labels
         cexCol=1.6,                      # magnification factor for column labels
         xlas=2,                          # column labels are to be perpendicular to the x-axis
         main="")                         # main titlek#legend(lx,0.97,legend=c("Not Transmitted","Transmitted"),bg="white",fill=transmissionTbl,title="transmission",cex=0.7)
#legend(lx,0.85,legend=c("Ng-Ct-","Ng-Ct+","Ng+Ct-","Ng+Ct+"),bg="white",fill=Ct_co_infectionTbl,title="Ct_co_infection",cex=0.7)
#legend(lx,0.38,legend=c("CEBA","permeability"),bg="white",fill=groupTbl,title="Group",cex=0.7)
legend(lx,0.97,legend=c("ATCC", "DESULPHO", "MD","PBS"),bg="white",fill=groupTbl,title="Group",cex=0.9)
legend(lx,0.75,legend=c("Colon","Jejunum"),bg="white",fill=siteTbl,title="Site",cex=0.9)
dev.off()

###~~~~~~~~set 3 fixed
mytotal2=read.csv("set3_mp_fixed.csv",header=T,row.names=1,check.names=F)

dim(mytotal2)
head(mytotal2)
mytotal[1:5,1:5]

mydata=mytotal
#mydata=t(mytotal)
mydata=mytotal[,c(-1:-3)] # remove first and second columns with all rows
mydata[1:5,1:5]

sitef=factor(mytotal$site)
levels(sitef)
siteTbl=c()
siteTbl[levels(sitef)[1]]="navy"
siteTbl[levels(sitef)[2]]="orange"

groupf=factor(mytotal$group)
levels(groupf)
groupTbl=c()
groupTbl[levels(groupf)[1]]="navy"
groupTbl[levels(groupf)[2]]="green"
groupTbl[levels(groupf)[3]]="red"
groupTbl[levels(groupf)[4]]="yellow"

#cs=colSums(JEdata)
#o=order(cs,decreasing=T)
#ct=JEdata[,o]
#ct[1:5,1:5]
#pt=t(apply(ct,1,function(x) 100*x/sum(x)))
pt=mydata #metaPhlAn result
pt[1:5,1:5]
nCols=25
dim(pt)
#hc.comp=hclust(dist(pt),method="comp")
#hc.ward=hclust(dist(pt),method="ward.D")

#nClrs=3
#memb.ward=cutree(hc.ward,k=nClrs)
#table(memb.ward)

lx=0.86
sideBars.ward=cbind(siteTbl[sitef],groupTbl[groupf])
colnames(sideBars.ward)=c("site","group")
#rownames(sideBars.ward)=memb.ward

#write.csv(file="cluster_ward.csv", memb.ward)

nCols=15
dim(pt)
postscript("heatmap_ward_CO_JE.eps",hor=T)

heatmap2(as.matrix(pt)[,1:nCols],         # show only nCols first columns of the table tbl
         col=rainbow(50,start=1/6,end=0),
         #col=greenred(100), # set heatmap color palette
         #col=colorRampPalette(c("darkblue","black","yellow"))(100),
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         Rowv = NA,
         #Rowv=as.dendrogram(hc.ward),  
         #Rowv=as.dendrogram(hc.ward),    # use our h. clustering
         RowSideColors=sideBars.ward,    # RowSideColors=sideBars,          # this sets three color columns: clustering, pH and Nugent score
         RowSideTitleCex=0.5,             # scaling factor for the titles of the color bars
         RowSideTitleLine=0.5,            # the distance of the titles for the color bars
         margins=c(18,4),                 # =c(bottom margin, right margin)
         #labRow=rownames(mydata),         # add row labels
         labRow=NA,                       # suppress row labels
         cexCol=0.8,                      # magnification factor for column labels
         xlas=2,                          # column labels are to be perpendicular to the x-axis
         main="")                         # main titlek#legend(lx,0.97,legend=c("Not Transmitted","Transmitted"),bg="white",fill=transmissionTbl,title="transmission",cex=0.7)
#legend(lx,0.85,legend=c("Ng-Ct-","Ng-Ct+","Ng+Ct-","Ng+Ct+"),bg="white",fill=Ct_co_infectionTbl,title="Ct_co_infection",cex=0.7)
#legend(lx,0.38,legend=c("CEBA","permeability"),bg="white",fill=groupTbl,title="Group",cex=0.7)
legend(lx,0.97,legend=c("ATCC", "DESULPHO", "MD","PBS"),bg="white",fill=groupTbl,title="Group",cex=0.9)
legend(lx,0.75,legend=c("Colon","Jejunum"),bg="white",fill=siteTbl,title="Site",cex=0.9)
dev.off()


