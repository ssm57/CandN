rm(list=ls(all=TRUE))   # clean environment
getwd()
#update.packages(ask = FALSE, dependencies = c('Suggests')) #update all installed R packages

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
library(scales)         # for points transparency on plots
library(gplots)         # for heatmap.2 function

setwd("/Users/Zireael/Desktop/Maslov/CN_paper") # replace with your working directory

# load some custom functions
source("src/functions.R")   

# read data 
states<-fread("data/2Cx2Nx4S_MC.txt", header = F, stringsAsFactors = F, data.table = F)
# separate flux and states tables
fluxes<-states[, 1:4]  
sub.states<-states[, 5:ncol(states)] 

##################################################
####            Prepare PCA (Fig 2A)          #### 
##################################################
# calc mean flux coord for each state
f.means<-data.frame()
for (i in 1:ncol(sub.states)){
  print(i)
  #i=1
  st.ids<-which(sub.states[,i]>0)
  sub.flux<-fluxes[st.ids,]
  f.means<-rbind(f.means, colMeans(sub.flux))
}

colnames(f.means)<-c("C1", "C2", "N1", "N2")

# compute PCA
pca<-prcomp(f.means)
summary(pca) # check - 1,2 PC -79% variance 

# plot PCA
tmp<-as.data.frame(pca$x[,1:2])             # get first two comp
tmp$names<-colnames(sub.states)             # if we want to display sample ids..

# calculate centroids for all states in C:N space
fc.means<-vector()
fn.means<-vector()
for (i in 1:ncol(sub.states)){
  print(i)
  #i=1
  st.ids<-which(sub.states[,i]==1)
  sub.flux<-fluxes[st.ids,]
  fc<-rowMeans(sub.flux[,1:2]) #1:3
  fn<-rowMeans(sub.flux[,3:4]) #4:6
  fc.means[i]<-mean(fc)
  fn.means[i]<-mean(fn)
}
rc.rn.rat<-fc.means/fn.means

# prepare colours palette by average fluxes
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(unique(rc.rn.rat)))
rats<-sort(unique(rc.rn.rat))
tmp$color<-"black"
for (i in 1:length(rats)){
  tmp$color[which(rc.rn.rat==rats[i])]<-my_palette[i]
}
# prepare symbols and colors for centroids
tmp$pch_col<-c("royal blue")
tmp$pch_col[27:33]<-"red3"
tmp$pch<-16
tmp$pch[33]<-4

par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)   # change graph params for legend 
xcoord<-jitter(tmp$PC1)#, amount = 5)      # preturb coords if needed
ycoord<-jitter(tmp$PC2)#, amount = 5)

cairo_pdf('graphs/pca_means_2x2_ell_25_corr.pdf', width = 7, height = 7) # save plot

# create plot grid
plot(x=xcoord, y=ycoord, col="white", pch=tmp$pch, cex=1, 
     xlab="PC1", ylab="PC2", 
     xlim=c(-600, 600), ylim=c(-600, 600))

### plot ellipses
library(car)
percent<-0.25
center.pch<-4
for (i in 1:ncol(sub.states)){ #i=1
  print(i)
  st.ids<-which(sub.states[,i]>0)
  sub.flux<-as.data.frame(fluxes[st.ids,])
  for (j in 1:ncol(sub.flux)){ # data centered
    sub.flux[,j]<-sub.flux[,j]-pca$center[j]
  }
  ttt<-as.matrix(sub.flux) %*% pca$rotation
  ttt.sub<-ttt[,1:2] 
  if (i<27){
    dataEllipse(ttt.sub[,1], ttt.sub[,2], levels=percent, 
                col = alpha(tmp$color[i], 0.2), add=T, plot.points=F,
                fill=T, fill.alpha=0.2, center.cex=0.1, center.pch=center.pch) #fill.alpha=0.1,
  }
  if (i>=27&i<33){
    dataEllipse(ttt.sub[,1], ttt.sub[,2], levels=percent, 
                col = alpha(tmp$color[i], 0.8), add=T, plot.points=F,
                fill=T, fill.alpha=0.2, center.cex=0.1, center.pch=center.pch)
    dataEllipse(ttt.sub[,1], ttt.sub[,2], levels=percent,
                col = alpha("black", 0.8), add=T, plot.points=F,
                fill=F, fill.alpha=0.1, center.cex=0.1, center.pch=center.pch)
  }
  if (i==33){
    dataEllipse(ttt.sub[,1], ttt.sub[,2], levels=percent, 
                col = alpha(tmp$color[i], 0.8), add=T, plot.points=F,
                fill=T, fill.alpha=0.2, center.cex=0.1, center.pch=center.pch, lty=2)
    dataEllipse(ttt.sub[,1], ttt.sub[,2], levels=percent,
                col = alpha("black", 0.8), add=T, plot.points=F,
                fill=F, fill.alpha=0.1, center.cex=0.1, center.pch=center.pch, lty=2)
  }
}

points(x=xcoord, y=ycoord, col=tmp$pch_col, pch=tmp$pch, cex=1) 

#### addition to print arrows for nutrient fluxes
arrow.len<-0.1                                # Length of the arrows about to plot.
choices<-1:2 
scores<-pca$x                                  # The scores
lam<-pca$sdev[choices]                        # Sqrt e-vals (lambda) 2 PC's
n<-nrow(scores)                               # no. rows scores
lam<-lam * sqrt(n)                            # See below.

y<-t(t(pca$rotation[,choices]) * lam)         # scaled eigenvecs (loadings)
# The scaled e-vecs are further reduced to 50% of their value
arrows(0, 0, y[, 1L] * 0.5, y[, 2L] * 0.5, length = arrow.len, col = "black")
ylabs<-dimnames(y)[[1L]]                      # Names of original coords
ylabs<-as.character(ylabs)
text(x = (y[, 1L]-60) * 0.5, y = (y[, 2L]-10) * 0.5, 
     labels = ylabs, col = "black", cex=0.8)   # Prints the names

# add labels for states
text(x=tmp$PC1[27:33], y=tmp$PC2[27:33], labels = c(1:length(tmp$names[27:33])), 
     pos = 3, cex=0.8, col="red3")
text(x=tmp$PC1[1:26], y=tmp$PC2[1:26], labels = c(7+1:length(tmp$names[1:26])), 
     pos = 3, cex=0.8, col="royal blue")

# legend for coloring according to <fc>/<fn>
legend.col(col = my_palette, lev = rats)
legend("topright", title="<f_c>/<f_n>", inset=c(-0.1,-0.07),
       bty="n", legend = c(""), lty= 0)
# legend for plot symbols
legend("topleft", inset=c(0,-0.15), bty="n", 
       legend = c("centroids of ISS", "centroids of UISS", "centroids of UIS"), 
       col = c("royal blue", "red3", "red3"), lty= 0, pch = c(16, 16, 4), cex=1)

dev.off() # end save plot

##################################################
####         Overlap network (Fig 3A)         ####
##################################################

# calculate state volumes
volumes<-colSums(sub.states)
uis.states<-sub.states[,27:33] # subset UIS

# calculate network of overlaps between UIS
overlaps<-matrix(0, ncol=ncol(uis.states), nrow=ncol(uis.states))
for (i in 1:nrow(overlaps)){
  for (j in i:ncol(overlaps)){
    overlaps[i,j]<-length(which(uis.states[,i]==1&uis.states[,j]==1))
    overlaps[j,i]<-overlaps[i,j]
  }
}

write.table(overlaps, "output/2Cx2Nx4S_UIS_overlaps.txt", quote=F, sep="\t",col.names = F, row.names = F)

##################################################
####         Yields and Multist (Fig 4A)      ####
##################################################
# load data for 4000 experiments
data<-fread("data_to_plot/Data_L2.txt", header = F, 
            stringsAsFactors = F, data.table = F)
# calc yield ratios
sds<-data.frame()
for (i in 1:nrow(data)){
  rat<-sd(log(data[i,1:4]/data[i,5:8]))
  sds<-rbind(sds, c(rat, data[i,11]/100000))
}

# create color palette and symbols for different distributions
my_palette <- colorRampPalette(rev(brewer.pal(9,"Greys")))(7)
#my_palette <- c('olivedrab3','firebrick3','royalblue','gold','magenta')
sds<-cbind(sds, c(rep(my_palette[3],1000), rep(my_palette[1],1000), 
                  rep(my_palette[2],1000), rep(my_palette[4],1000)))

sds<-cbind(sds, c(rep(16,1000), rep(0,1000), 
                  rep(17,1000), rep(18,1000)))

