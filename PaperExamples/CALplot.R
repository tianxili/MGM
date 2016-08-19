library(igraph)

X <- read.table("CAL500_stability_result2.txt",sep=",")

X <- as.matrix(X)
n <- nrow(X)
A <- matrix(0,n,n)
A[X>99] <- 1
sum(A)
g <- graph.adjacency(A,mode="undirected")
V(g)$shape = c(rep("circle",118),rep("square",16))
V(g)$label.cex = 0.5

node.class <- rep("",118)
node.class[1:36] <- "emotions"
node.class[37:55] <- "genre"
node.class[56:73] <- "instrument"
node.class[74:98] <-"song"
node.class[99:109] <- "usage"
node.class[110:118] <- "vocals"
node.class[119:122] <- "mean_timbre"
node.class[123:126] <- "std_timbre"
node.class[127:130] <- "mean_std_timbre"
node.class[131:134] <- "std_std_timbre"
node.color.factor <- factor(node.class)
r.color <- rainbow(15)
r.color <- r.color[-c(9,10)]
V(g)$color <- r.color[
as.numeric(node.color.factor)]
legend.color <- V(g)$color[c(1,37,56,74,99,110,119,123,127,131)]
legend.text <- node.class[c(1,37,56,74,99,110,119,123,127,131)]
legend.shape <- c(rep(21,6),rep(22,4))
V(g)$label <- c(1:118,1:16)
lay <- layout.fruchterman.reingold(g)

plot(g,vertex.size=8,layout=lay,margin=c(0,0.5,0,0.5))
legend("topleft",legend.text,pch=legend.shape,col=legend.color)
