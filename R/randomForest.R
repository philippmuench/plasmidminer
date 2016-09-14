library(ranger)

pos <- as.data.frame(read.csv("../pos.csv", header=T,sep=";"))
pos$name <- NULL
pos$id <- NULL
pos$length <- NULL
# unfactorize
pos[,colnames(pos)]<-lapply(colnames(pos), function(x) as.numeric(as.character(pos[,x])))
pos$type <- "pos"

neg <- as.data.frame(read.csv("../neg.csv", header=T,sep=";"))
neg$name <- NULL
neg$id <- NULL
neg$length <- NULL
# unfactorize
neg[,colnames(neg)]<-lapply(colnames(neg), function(x) as.numeric(as.character(neg[,x])))
neg$type <- "neg"

neg <- neg[1:nrow(pos),]

training <- rbind(pos, neg)
training <- training[complete.cases(training),]
training$type <- as.factor(as.character(as.matrix(training$type)))
forest <- ranger(type ~ ., data = training, num.trees=1000, num.threads=2)


library(rpart)
library(rpart.plot)
tree <- rpart(type ~ .,data=training)

# plot tree
rpart.plot(tree)