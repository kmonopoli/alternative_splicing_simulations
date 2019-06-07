# install.packages("XQuartz") # need to run this first to install summarytools
# install.packages( pkgs = "summarytools" )
library(summarytools)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(reshape)
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


# Counts
ct_dat <- read.csv(file=paste0(DIR,"export_start_posns.csv"), header=TRUE, sep=",")
Molten <- melt(ct_dat, id.vars = "index")
hplt2 <- ggplot(Molten, aes(x = index, y = value, colour = variable))+
      # geom_line()+
      geom_point()+xlab("Start Position")+ylab("Count")

#OLD histogram
# h_dat <- read.csv(file=paste0(DIR,"export_hist_starts.csv"), header=TRUE, sep=",")
# Molten <- melt(h_dat, id.vars = "bin_start")
# hplt2 <- ggplot(Molten, aes(x = bin_start, y = value, colour = variable))+ 
#       geom_line()+geom_point()+xlab("Start Position")+ylab("Frequency")



# Full Transcripts
# get other dataframe output
ft <- read.csv(file=paste0(DIR,"export_full_transcripts.csv"), header=TRUE, sep=",")
ft <-ft[order(ft$splice_type_read),] 
x_1 =ft$start_1 
xend_1 = ft$end_1
x_2 =ft$start_2 
xend_2 = ft$end_2
x_3 =ft$start_3
xend_3 = ft$end_3

x = c(x_1,x_2,x_3)
xend = c(xend_1,xend_2,xend_3)
y =c(seq(1, length(x_1), by=1),seq(1, length(x_1), by=1),seq(1, length(x_1), by=1))
spltype= c(as.character(ft$splice_type_read),as.character(ft$splice_type_read),as.character(ft$splice_type_read))

segment_data = data.frame(
  x = x,
  xend = xend,
  y =y,
  yend = y,
  splice_type = spltype
)

txplt <- ggplot(segment_data, aes(x = x, y = y, xend = xend, yend = yend,colour=splice_type))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = y))

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

y_height = max(segment_data$y)+max(segment_data$y)/10
gy = c(y_height,y_height,y_height,y_height,y_height)


segment_data2 = data.frame(
  x = gx,
  xend = gxend,
  y =gy,
  yend = gy,
  splice_type = c("exon","intron","exon","intron","exon")
)
sz <- c(5,1,5,1,5)

txplt2 <- ggplot(segment_data2, aes(x = gx, y = gy, xend = gxend, yend = gy))+
  geom_segment(aes(x = gx, y = gy, xend = gxend, yend = gy,size=sz))

txplt_all <- ggplot(NULL, aes(x = x, y = y, xend = xend, yend = yend,colour=splice_type)) + 
    geom_segment(data = segment_data) +
    geom_segment(data = segment_data2, size = sz, color = c("black","black","black","black","black"))+
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
    labs("exon3 included","exon3 excluded","unspliced","unspliced (short)")


# Plot chart and table into one object
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             c(3,3,3))

grid.arrange(txplt_all, hplt2+theme(legend.position="none"), tbl, layout_matrix = lay)







