library(lattice)

data <- read.table("data.txt", sep="\t", header=F, stringsAsFactors=F)

pdf("rarefaction.pdf");
xyplot(data$V3 ~ data$V2, group=data$V1, data, type="l", xlim=c(0,max(data$V2)), ylim=c(0, max(data$V3)), xlab="Number of Reads", ylab="Species Count", auto.key=T); # par.settings=list(superpose.line=list(col=c("red","blue"))));

dev.off();
