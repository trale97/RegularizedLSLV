#### Title: Empirical Data Application
#### Author: Tra Le 
#### Supervisor: Dr. Katrijn Van Deun
#### Created: May 30, 2023
#### Last modified: June 10, 2023

#########################################################################################
#####                             GENE EXPRESSION                                   #####
#########################################################################################

rm(list = ls(all.names = TRUE))

### Download data set 
Data_loacation= 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE7nnn/GSE7329/matrix/GSE7329_series_matrix.txt.gz'
destfile= paste0(getwd(),"/GSE7329_series_matrix.txt.gz")
download.file(Data_loacation,destfile)

full_data=gzfile('GSE7329_series_matrix.txt.gz','rt')  
full_data=read.csv(full_data,header=F,sep = '\t')

col_names= c()
for (i in 736:751) {
  # col_names=c(col_names, as.character(full_data[i,1]),as.character(dat[i,2]))
  col_names=c(col_names,sapply(full_data[i,], as.character))
} 
View(col_names)

character <- c()
for (j in 287:302){
  character <- c(character,sapply(full_data[j,], as.character))
}
View(character)
length(character)
character <- character[-c(1,32)]
character <- character[-c(12, 13, 27)]
character_new <- c(rep("dup15",7), rep("FMR1", 6), rep("control", 14))

matrix1= full_data[752:703647,]
matrix2= matrix(rep(0, prod(dim(matrix1))), ncol = dim(matrix1)[2])
for(i in 1:nrow(matrix1) ){
  matrix2[i,]= as.numeric(sapply(matrix1[i,],as.character))
}

Data_Autism=matrix(data = t(matrix2),byrow = TRUE, ncol = length(col_names))
colnames(Data_Autism)=col_names # adding names
IndexCol= Data_Autism[,1]
Data_Autism=Data_Autism[,-c(1,dim(Data_Autism)[2])] # removing first and last column 
Data_Autism=Data_Autism[,-c(12,13,27)] # Removing mislabel individuals 
dim(Data_Autism)

### Removing NA values 
Data_Autism= t(Data_Autism)
Data_Autism=scale(Data_Autism, center = TRUE, scale = TRUE)

NA_Data_Autism= is.na(Data_Autism)
NA_index= colSums(NA_Data_Autism)!=0
sum(NA_index) # number of NA columns
which(NA_index!=0)
Data_Autism=Data_Autism[,-which(NA_index!=0)]
dim(Data_Autism)

Data_Autism=scale(Data_Autism, center = TRUE, scale = TRUE)
save(Data_Autism, file = "Data_Autism_scaled.RData")
load("Data_Autism_scaled.RData")

### Analyze using CCLSLV
X <- Data_Autism
J= dim(X)[2]
cardinality= floor(seq(from=100,to=(J-100),length.out=100))*3
IS_genetic <- matrix(ncol = 3, nrow = length(cardinality))
for (i in length(cardinality)){
  a <- IS_CCLSLV(X, R = 3, card = cardinality[i], MaxIter = 200, eps = 10^-6, nstarts = 5)
  IS_genetic[i,] <- c(cardinality[i], a$value, a$vaf)
}
genetic_model <- MULTISTART_CCLSLV(X, R = 3, CARD = cardinality[which.max(IS_genetic[,2])], MaxIter = 200, eps = 10^-6, nstarts = 10)
componentscores <- genetic_model$scores
componentscores <- cbind(componentscores, character_new)
componentscores <- as.data.frame(componentscores)
colnames(componentscores) <- c("LV_1", "LV_2", "LV_3", "Group")
componentscores <- as.data.frame(componentscores)
### 3D plot
#install.packages("plotly")
library(plotly)

p <- plot_ly(componentscores, x =~ Component_1, y =~ Component_2, z =~ Component_3,
             color =~ Group, colors = c("red", "blue", "#32a852")) %>%
  add_markers(size = 1)

### Analyze using RLSLV (have to try several times with different lambda)
RLSLV_genetic <- MULTISTART(X, R = 3, MaxIter = 200, eps = 10^-6, nstarts = 10, lambda = 16.85)
sum(RLSLV_genetic$loadings != 0)
dim(RLSLV_genetic$loadings)
