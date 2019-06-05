library(ggplot2)
# Get data
DIR = "/Users/kmonopoli/Google_Drive/research/alt_splicing_simulations/as_simulations_km/"
dat <- read.csv(file=paste0(DIR,"export_dataframe.csv"), header=TRUE, sep=",")
names(dat)[names(dat) == "start"] <- "seq_start"

# Scatter plot with the 2d density estimation
sp <- ggplot(dat, aes(x=splice_type, y=seq_start))  + geom_bin2d()+ scale_x_continuous(breaks=c(0,1,2,3,4),labels=c("unspliced", "spliced", "alt spliced","alt spliced (short)","unspliced (short)"), limits = c(0,4.5))

# Rename / Change tick mark labels
# Create summary table
SummaryTable <- data.frame(
  h_ae = dat$h_ae[1],
  h = dat$h[1],
  exon_se =dat$exon_se[1],
  exon = dat$exon[1],
  trnsc_rate = dat$transc_rate[1],
  intron = dat$intron[1],
  intron_1 = dat$intron_1[1],
  intron_2 = dat$intron_2[1]
)

# Create a table plot
library(gridExtra)
mytheme <- ttheme_default(base_size = 12, colhead=list(fg_params = list(parse=TRUE)), padding = unit(c(5, 5), "mm") )
tbl <- tableGrob(SummaryTable, rows=NULL, theme = mytheme)



# histogram
hplt<-ggplot(dat, aes(x=seq_start, color=splice_type_read, fill=splice_type_read)) + 
  geom_histogram( position="dodge")# +

# junction reads
# get only junction read data
jrds <-subset(dat, dat$junction != "no")
jplt<- ggplot(jrds, aes(x=seq_start, color=junction, fill=junction)) + 
  geom_histogram( position="dodge")

# Plot chart and table into one object
lay <- rbind(c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(2,2,2,3,3),
             c(2,2,2,3,3),
             c(4,4,4,4,4))

grid.arrange(sp, hplt, jplt, tbl, layout_matrix = lay)
