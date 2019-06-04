
library(ggplot2)

DIR = "/Users/kmonopoli/Google_Drive/research/alt_splicing_simulations/as_simulations_km/"
dat <- read.csv(file=paste0(DIR,"export_dataframe.csv"), header=TRUE, sep=",")
names(dat)[names(dat) == "start"] <- "seq_start"
#ggplot(dat, aes(x=splice_type, y=seq_start)) + geom_point()


# Scatter plot with the 2d density estimation
sp <- ggplot(dat, aes(x=splice_type, y=seq_start)) #+geom_point()
sp + geom_bin2d()
