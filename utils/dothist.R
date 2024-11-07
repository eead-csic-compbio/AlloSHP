
library(stringr)

args <- commandArgs(TRUE)

dotfile1 <- readLines(args[1]) #"Bdis.fna.gz.Bsta.fna.gz_GSAlign_-no_vcf_-one.dot")
dotfile2 <- readLines(args[2]) #"Bdis.fna.gz.Bsta.fna.gz_Cgaln_-K11_-BS10000_-X12000_-fc_-cons.dot")

len1 <- str_match(dotfile1, 'length: (\\d+)')
hist1 <- hist(as.numeric(len1[,2]))

len2 <- str_match(dotfile2, 'length: (\\d+)')
hist2 <- hist(as.numeric(len2[,2]))
maxf <- max(hist1$counts,hist2$counts)

col1 = rgb(0,0,1,1/4)
col2 = rgb(1,0,0,1/4)

plot( hist1, col=col1, ylim=c(0,maxf), main="", xlab="Alignment length", )  
plot( hist2, col=col2, ylim=c(0,maxf), add=T) 

legend("topright", c(args[1], args[2]), col=c(col1, col2), lwd=10)
