# install.packages("XQuartz") # need to run this first to install summarytools
# install.packages( pkgs = "summarytools" )
# install.packages("wesanderson")
# install.packages("hexbin")
library(summarytools)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(reshape)
library("wesanderson")
library(RColorBrewer)
library(grid)
library(stringr)
# library(hexbin)

# Get data
DIR = "/Users/kmonopoli/Google_Drive/research/alt_splicing_simulations/as_simulations_km/"
dat <- read.csv(file=paste0(DIR,"sim_outputs/export_dataframe_1.csv"), header=TRUE, sep=",")
names(dat)[names(dat) == "start"] <- "seq_start"

tbl_nms = c("half life intron 1 (min)","half life intron 2 (min)","half life intron 3 (min)",
            "u dist (nt)","alt exon len (nt)","exon len (nt)",
            "transc rate (nt/min)","intron 1 len (nt)",
            "intron 2 len (nt)", "intron 3 len (nt)","psi alt exon")
# Create summary table
SummaryTable <- data.frame(
  h_intron_1 = dat$h_intron_1[1],
  h_intron_2 = dat$h_intron_2[1],
  h_intron_3 = dat$h_intron_3[1],
  u_dist = dat$u_dist[1],
  exon_se =dat$exon_se[1],
  exon = dat$exon[1],
  trnsc_rate = dat$transc_rate[1],
  intron_1 = dat$intron_1[1],
  intron_2 = dat$intron_2[1],
  intron_3 = dat$intron_3[1],
  psi_se = dat$psi_se[1]
)
library(gridExtra)
mytheme <- ttheme_default(base_size = 8, colhead=list(fg_params = list(parse=TRUE)), padding = unit(c(5, 5), "mm") )
tbl <- tableGrob(SummaryTable, cols=tbl_nms,rows=NULL, theme = mytheme)



# FREQUENCY DIAGRAM
# counts
ct_dat <- read.csv(file=paste0(DIR,"sim_outputs/export_start_posns_1.csv"), header=TRUE, sep=",")
max_x <- round(max(ct_dat$index)/100,digits=0)*100
min_x <- round(min(ct_dat$index)/100,digits=0)*100

ct_dat<- ct_dat[colSums(!is.na(ct_dat)) > 0] # drop any columns that are empty (all NA's)
names(ct_dat)[names(ct_dat) == "ct_intron_1_excluded"] <- "intron 1 excluded"
names(ct_dat)[names(ct_dat) == "ct_intron_1_intron_2_excluded"] <- "intron 1 & 2 excluded"
names(ct_dat)[names(ct_dat) == "ct_intron_2_excluded"] <- "intron 2 excluded"
names(ct_dat)[names(ct_dat) == "ct_intron_3_excluded"] <- "intron 3 excluded"
names(ct_dat)[names(ct_dat) == "ct_unspliced"] <- "unspliced"

Molten <- melt(ct_dat, id.vars = "index")
Molten$variable <- factor(Molten$variable, levels(Molten$variable)[c(2,1,3:5)])
# Molten$variable <- factor(Molten$variable, levels(Molten$variable)[c(1,2,3:5)]) # if only 2 vars
# sort by variable (flip so unspliced at bottom)
hplt <- ggplot(Molten, aes(x = index, y = value, colour = (Molten$variable)))+
      geom_point(alpha = 0.5)+#,shape = ".")+
      #geom_line()+geom_path()+
      guides( colour = guide_legend(override.aes = list(alpha=1)))+
      theme(legend.position="bottom")+
      xlim(min_x-1, max_x)+
      # scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
      xlab("Start Position")+ylab("Count")#+scale_color_manual(values=brewer.pal(n=ncol(ct_dat)-1, name="Set2"))

# GENE DIAGRAM
# Full Transcripts
# get other dataframe output
ft <- read.csv(file=paste0(DIR,"sim_outputs/export_full_transcripts_1.csv"), header=TRUE, sep=",")

ft <-ft[order(ft$splice_type_read,ft$lengths),]

x_1 =ft$start_1-ft$u_dist[1]
xend_1 = ft$end_1-ft$u_dist[1]
x_2 =ft$start_2-ft$u_dist[1]
xend_2 = ft$end_2-ft$u_dist[1]
x_3 =ft$start_3-ft$u_dist[1]
xend_3 = ft$end_3-ft$u_dist[1]

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
e_2 = iend_1
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
  x = gx-ft$u_dist[1],
  xend = gxend-ft$u_dist[1],
  y =gy,
  yend = gy,
  splice_type = c("exon","intron","exon","intron","exon")
)
nms <- c("u dist","intron 1","exon se", "intron 2", "exon")
sz <- c(5,1,5,1,5)

# txplt2 <- ggplot(segment_data2, aes(x = x, y = y, xend = xend, yend = y))+
#   geom_segment(aes(x = x, y = y, xend = xend, yend = y,size=sz))
#   geom_text(data = segment_data2, aes(label=splice_type), position=position_nudge(x=0,y=1), hjust = -0.5)


txplt_all <- ggplot(NULL, aes(x = x, y = y, xend = xend, yend = yend,colour=splice_type)) +
    geom_segment(data = segment_data) +
    geom_segment(data = segment_data2, size = sz, color = c("black","black","black","black","black"))+
    # theme(legend.position="bottom")+#, axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
    # scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0, 0))+
    xlim(min_x-1, max_x)+
    geom_text(data = segment_data2, aes(label=nms), position=position_nudge(x=0,y=600), hjust = 0, size = 2.8, colour = "black")
  




txplt_all <- txplt_all+theme(legend.position="none", axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())
  
hplt <-hplt +theme(legend.title = element_blank(),legend.position="bottom",legend.spacing.x = unit(1.0, 'cm'))

gA <- ggplotGrob(txplt_all)
gB <- ggplotGrob(hplt)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
# grid.arrange(gA, gB, ncol=1)




## Junction Reads
# exon-exon junction reads
jnc_dat <- read.csv(file=paste0(DIR,"sim_outputs/export_junction_reads_1.csv"), header=TRUE, sep=",")
jnc_dat <-jnc_dat[order(jnc_dat$junction.read,jnc_dat$start.position),]
jnc_dat<-jnc_dat[!(jnc_dat$junction.read=="not junction read"),] # remove no's

x_1 =jnc_dat$start.position
xend_1 = jnc_dat$start.position + jnc_dat$read.length

y_1 =c(seq(1, length(x_1), by=1))
junct_read= c(as.character(jnc_dat$junction.read))
junct_read_type= c(as.character(jnc_dat$junction_read_type))

junc_segment_data = data.frame(
  x = x_1,
  xend = xend_1,
  y =y_1,
  yend = y_1,
  jnct_read = junct_read,
  jnct_read_type = junct_read_type
)

  
jncplt <- ggplot(junc_segment_data, aes(x = x, y = y, xend = xend, yend = yend,colour = jnct_read_type))+#colour=jnct_read))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = y))




y_height_jnc= max(junc_segment_data$y)/2
segment_data3 <- segment_data2
segment_data3$y <-c(y_height_jnc,y_height_jnc,y_height_jnc,y_height_jnc,y_height_jnc)
segment_data3$yend <-c(y_height_jnc,y_height_jnc,y_height_jnc,y_height_jnc,y_height_jnc)
# segment_data3$x <- segment_data3$x - dat$u_dist[1]
# segment_data3$xend <- segment_data3$xend - dat$u_dist[1]
sz2 <- c(300,1,300,1,300)

txplt_junc <- ggplot(NULL, aes(x = x, y = y, xend = xend, yend = yend,colour = jnct_read))+#colour=jnct_read)) +
  geom_segment(data = segment_data3, size = sz2, color = c("grey","grey","grey","grey","grey"))+
  geom_text(data = segment_data3, aes(label=nms), position=position_nudge(x=0,y=300), hjust = 0, size = 2.8, colour = "black")+
  geom_segment(data = junc_segment_data) +
  xlab("Start Position") +
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.title = element_blank())
  





# Plot chart and table into one object
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2),
             # c(3,3,3),
             # c(3,3,3),
             # c(3,3,3),
             c(4,4,4))

# grid.arrange(gA, gB,
#              # txplt_junc,
#              tbl, layout_matrix=lay)


# Plot junction reads and table into one object
lay2 <- rbind(c(1,1,1),
             c(1,1,1),
             c(2,2,2))

grid.arrange(txplt_junc,
             tbl, layout_matrix=lay2)
