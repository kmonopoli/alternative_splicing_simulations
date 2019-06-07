# install.packages("XQuartz") # need to run this first to install summarytools
# install.packages( pkgs = "summarytools" )
library(summarytools)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
# Get data
DIR = "/Users/kmonopoli/Google_Drive/research/alt_splicing_simulations/as_simulations_km/"
dat <- read.csv(file=paste0(DIR,"export_dataframe.csv"), header=TRUE, sep=",")
names(dat)[names(dat) == "start"] <- "seq_start"


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
library(gridExtra)
mytheme <- ttheme_default(base_size = 12, colhead=list(fg_params = list(parse=TRUE)), padding = unit(c(5, 5), "mm") )
tbl <- tableGrob(SummaryTable, rows=NULL, theme = mytheme)



# histogram
hplt<-ggplot(dat, aes(x=read_start, color=splice_type_read, fill=splice_type_read)) + 
  geom_histogram( position="dodge")







# TODO:
# Scatter plot with the 2d density estimation
# sp <- ggplot(dat, aes(x=splice_type, y=seq_start))  + geom_bin2d()+ scale_x_continuous(breaks=c(0,1,2,3,4),labels=c("unspliced", "spliced", "alt spliced","alt spliced (short)","unspliced (short)"), limits = c(0,4.5))

# Rename / Change tick mark labels






# Full Transcripts
# get other dataframe output

ft <- read.csv(file=paste0(DIR,"export_full_transcripts.csv"), header=TRUE, sep=",")

x_1 =ft$start_1 #c(1, 1)
xend_1 = ft$end_1#c(5,10)
x_2 =ft$start_2 
xend_2 = ft$end_2
x_3 =ft$start_3
xend_3 = ft$end_3





#x = c(gx,x_1,x_2,x_3)
x = c(x_1,x_2,x_3)

# xend = c(gxend,xend_1,xend_2,xend_3)
xend = c(xend_1,xend_2,xend_3)

# y =c(gy,seq(1, length(x_1), by=1),seq(1, length(x_1), by=1),seq(1, length(x_1), by=1))
y =c(seq(1, length(x_1), by=1),seq(1, length(x_1), by=1),seq(1, length(x_1), by=1))

# spltype= c(rep("gene model",length(gx)),as.character(ft$splice_type_read),as.character(ft$splice_type_read),as.character(ft$splice_type_read))
spltype= c(as.character(ft$splice_type_read),as.character(ft$splice_type_read),as.character(ft$splice_type_read))

segment_data = data.frame(
  x = x,
  xend = xend,
  y =y,
  yend = y,
  splice_type = spltype
)

txplt <- ggplot(segment_data, aes(x = x, y = y, xend = xend, yend = yend,colour=splice_type))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))#+
  #geom_segment(aes(x = gx, y = gy, xend = gxend, yend = gy))

# gene diagram
e_1 = c(0)
eend_1 = c(ft$u_dist[1])
i_1 = eend_1
iend_1 =  c(ft$u_dist[1]+ft$intron_1[1])
e_2 = i_end
eend_2 = c(ft$u_dist[1]+ft$intron_1[1]+ft$exon_se[1])
i_2 = eend_2
iend_2 = c(ft$u_dist[1]+ft$intron_1[1]+ft$exon_se[1]+ft$intron_2[1])
e_3 = iend_2
eend_3 = c(ft$u_dist[1]+ft$intron_1[1]+ft$exon_se[1]+ft$intron_2[1]+ft$exon[1])

gx = c(e_1,i_1,e_2,i_2,e_3)
gxend=c(eend_1,iend_1,eend_2,iend_2,eend_3)
gy = c(500,500,500,500,500)



segment_data2 = data.frame(
  x = gx,
  xend = gxend,
  y =gy,
  yend = gy,
  splice_type = c("exon","intron","exon","intron","exon")
)
sz <- c(5,2,5,2,5)

txplt2 <- ggplot(segment_data2, aes(x = gx, y = gy, xend = gxend, yend = gy))+
  geom_segment(aes(x = gx, y = gy, xend = gxend, yend = gy,size=sz))#+

txplt_all <- ggplot(NULL, aes(x = x, y = y, xend = xend, yend = yend,colour=splice_type)) + 
    geom_segment(data = segment_data) +
    geom_segment(data = segment_data2, size = sz, color = c("black","black","black","black","black"))



# Plot chart and table into one object
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             c(3,3,3))

grid.arrange(txplt_all, hplt, tbl, layout_matrix = lay)






















# 
# # density plot
# dplt<-ggplot(dat, aes(x=read_start, color=splice_type_read)) + 
#   geom_density()+scale_fill_brewer(palette="Paired")
# 
# dplt
# 
# 
# 
# 
# 
# # dot histogram
# # get frequencies for each type of splicing (splice_type_read) for each start
# #   computes frequences (which outputs a matrix)
# #   deletes columns 2-5 of matrices
# #   turns matrix to dataframe and removes Total and <NA>
# freq_nothing_spliced <- head(melt(freq(dat[dat$splice_type_read == 'unspliced',]$read_start)[,-c(2,3,4,5)]),-2)
# freq_intron_3_excluded <- head(melt(freq(dat[dat$splice_type_read == 'spliced',]$read_start)[,-c(2,3,4,5)]),-2)
# freq_intron_1_2_excluded<- head(melt(freq(dat[dat$splice_type_read == 'alt spliced',]$read_start)[,-c(2,3,4,5)]),-2)
# freq_intron_1_excluded<- head(melt(freq(dat[dat$splice_type_read == 'alt spliced (short)',]$read_start)[,-c(2,3,4,5)]),-2)
# freq_nothing_spliced_incomplete <- head(melt(freq(dat[dat$splice_type_read == 'unspliced (short)',]$read_start)[,-c(2,3,4,5)]),-2) 
# 
# hplt1<-ggplot(freq_nothing_spliced_incomplete, aes(x= as.numeric(row.names(freq_nothing_spliced_incomplete)),y=value)) + 
#   geom_point(size=2, shape=23)
# hplt1
# 
# 
# # Test with freq_intron_3_excluded
# br = seq(min(dat$read_start),max(dat$read_start)+1000,by=100)
# freq_3 <-hist(dat[dat$splice_type_read == 'spliced',]$read_start, breaks=seq(0,130,by=1))
# freq_3
# 
# 
# # Plot chart and table into one object
# lay <- rbind(c(1,1),
#              c(1,1),
#              c(2,2),
#              c(2,2),
#              c(3,3))
# 
# grid.arrange(hplt, hplt1, tbl, layout_matrix = lay)
# 
# 
# 
# # # junction reads
# # # get only junction_ee read data
# # jeerds <-subset(dat, dat$junction_ee != "no")
# # jeeplt<- ggplot(jrds, aes(x=read_start, color=junction_ee, fill=junction_ee)) + 
# #   geom_histogram( position="dodge")
# # # get only junction_ie read data
# # jierds <-subset(dat, dat$junction_ie != "no")
# # jieplt<- ggplot(jrds, aes(x=read_start, color=junction_ie, fill=junction_ie)) + 
# #   geom_histogram( position="dodge")
# 
# 
# 
#
# gy = c(seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1),
#        seq((length(ft$exon_se)+150), (length(ft$exon_se)+250), by=1),
#        seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1),
#        seq((length(ft$exon_se)+150), (length(ft$exon_se)+250), by=1),
#        seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1))
# 
# gx=c(seq(e_1,length(seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1))),
#      seq(i_1,length(seq((length(ft$exon_se)+150), (length(ft$exon_se)+250), by=1))),
#      seq(e_2,length(seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1))),
#      seq(i_2,length(seq((length(ft$exon_se)+150), (length(ft$exon_se)+250), by=1))),
#      seq(e_3,length(seq((length(ft$exon_se)+100), (length(ft$exon_se)+300), by=1))))
