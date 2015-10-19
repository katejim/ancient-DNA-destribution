# load files in folder and return list of loaded files
loadGroupFiles <- function(path){
  files <- list.files(path, full.names = TRUE)
  filelist <- lapply(files, read.table, stringsAsFactors = FALSE, header = TRUE)
  filelist
}

# filter and combine data
readGroup <- function(in.files, aval.positions){
  files <- in.files
  result <- files[[1]]
  result <- subset(result, result$POS %in% aval.positions)
  result <- result[, 4:ncol(result)]
  
  i<-1
  for (i in 2:length(files)){
    temp <- files[[i]]
    temp <- subset(temp, temp$POS %in% aval.positions)
    temp <- temp[, 4:ncol(temp)]
    result <- cbind(result, temp)
    i <- i + 1
  }
  # filter by chromosome 1 position
  result
}



AFR <- loadGroupFiles("../data/VEC_1KG_v2/AFR")
EUR <- loadGroupFiles("../data/VEC_1KG_v2/EUR")
EAS <- loadGroupFiles("../data/VEC_1KG_v2/EAS")
AMR <- loadGroupFiles("../data/VEC_1KG_v2/AMR")
SAS <- loadGroupFiles("../data/VEC_1KG_v2/SAS")



#get unique elements in data column 
unique.elements <- function(data, column){
  rapply(unique(data[column]), c)
}

#use for list of files from directory, return common values in current column
common.values <- function(dataVector, column) {
  uniqueListChroms <- lapply(dataVector, unique.elements, column)
  Reduce(intersect, uniqueListChroms)
}

#use vector of list of files from directory, return common values for all files
common.modern <- function(data, column){
  intersect.in.pop <- lapply(data, common.values, column)
  Reduce(intersect, intersect.in.pop)
}

ancient <- loadGroupFiles("../data/ancient")
modern <- list(AFR, EUR, EAS, AMR, SAS)

#get common positions in modern
common.modern.pos <- common.modern(modern, 2)
#get common positions in ancient
common.ancient.pos <- common.values(ancient, 3)
#intersect
common.poss <- intersect(common.modern.pos, common.ancient.pos)


#filter populations
AMRR <- readGroup(AMR, common.poss)
AFRR <- readGroup(AFR, common.poss)
EURR <- readGroup(EUR, common.poss)
EASR <- readGroup(EAS, common.poss)
SASR <- readGroup(SAS, common.poss)

#filter ancient
Neanderthal <- ancient[[1]]
Neanderthal <- subset(Neanderthal, Neanderthal[[3]] %in% common.poss)
Chimpanzee <- ancient[[2]]
Chimpanzee <- subset(Chimpanzee, Chimpanzee[[3]] %in% common.poss)
Denisova <- ancient[[3]]
Denisova <- subset(Denisova, Denisova[[3]] %in% common.poss)

# Preparing for PCA Neanderthal + Denisova Chimpanzee
data <- cbind(Neanderthal[,4:ncol(Neanderthal)], Chimpanzee[,4:ncol(Chimpanzee)], Denisova[,4:ncol(Denisova)])
data <- data[seq(1,nrow(data), by = 1),]

# PCA Neanderthal + Denisova + Chimpanzee
pca <- prcomp(t(data), retx = TRUE)

#project modern people to ancient pca components
pred.AMR <- predict(pca, t(AMRR))
pred.AFR <- predict(pca, t(AFRR))
pred.EUR <- predict(pca, t(EURR))
pred.EAS <- predict(pca, t(EASR))
pred.SAS <- predict(pca, t(SASR))


colours <- c(1, 2, 3, 4, 5)
pchs <- c(15:19)

#plot ancient prople with modern projection
plot(pca$x[,1], pca$x[,2], col=colours, cex=1, xlab = paste0('PC1 ', " (", round(pca$sdev[1]/sum(pca$sdev)*100,0), '%)'),
ylab = paste0('PC2 ', " (", round(pca$sdev[2]/sum(pca$sdev)*100,0), '%)'),pch=19)

legend(x="topright", legend=c("Neanderthal", "Denisova", "Chimpanzee"), 
       col=colours,lty = c(2, -1, 1), pch = pchs,
       merge = TRUE, bg = "gray90")


points(pred.AFR[,1], pred.AFR[,2], col=colours[1], pch=pchs[1])
points(pred.EUR[,1], pred.EUR[,2], col=colours[2], pch=pchs[3])
points(pred.EAS[,1], pred.EAS[,2], col=colours[4], pch=pchs[4])
points(pred.AMR[,1], pred.AMR[,2], col=colours[3], pch=pchs[2])
points(pred.SAS[,1], pred.SAS[,2], col=colours[5], pch=pchs[5])




#plot only modern projection
plot(pred.AFR[,1], pred.AFR[,2], col=colours[1], pch=pchs[1])
points(pred.EUR[,1], pred.EUR[,2], col=colours[2], pch=pchs[2])
points(pred.EAS[,1], pred.EAS[,2], col=colours[4], pch=pchs[3])
points(pred.AMR[,1], pred.AMR[,2], col=colours[3], pch=pchs[4])
points(pred.SAS[,1], pred.SAS[,2], col=colours[5], pch=pchs[5])


legend(x="topright", legend=c("AFR",  "EUR",  "EAS", "AMR", "SAS"), 
       col=colours,lty = c(2, -1, 1), pch = pchs,
       merge = TRUE, bg = "gray90")
  


