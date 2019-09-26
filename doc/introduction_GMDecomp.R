## ---- eval = TRUE--------------------------------------------------------
library(GMDecomp)
load("tobacco_clr.Rdata")

## ---- eval = TRUE--------------------------------------------------------
tobacco.gmd <- GMD(X = tobacco_clr$data, H = tobacco_clr$H, Q = diag(1, dim(tobacco_clr$data)[2]), K = 10)

## ---- eval = TRUE, fig.height= 6, fig.width= 6---------------------------
screeplot(tobacco.gmd) #the screeplot of the top 10 GMD components 

## ---- eval = TRUE, fig.height= 6, fig.width= 6---------------------------
gmd.order = order(rowSums(tobacco.gmd$V[,1:2]^2), decreasing = T)
plot.index = gmd.order[1:3]
plot.names = tobacco_clr$otu.names[plot.index]
biplot(fit = tobacco.gmd, index = plot.index, names = plot.names, sample.col = tobacco_clr$sample.color, sample.pch = tobacco_clr$sample.pch, arrow.col = 'grey50')

