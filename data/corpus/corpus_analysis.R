# Corpus analysis for Cogsci2014 Referential Complexity paper (Lewis & Frank)
#
# This script reproduces the statistical analyses for all four corpora (part 1), and 
# reproduces the plot presented in the paper (part 2).

rm(list=ls())

library(RColorBrewer)

# define functions
bm.partial<-function(x,y,z) {round((cor(x,y, use="complete.obs")-cor(x,z, use="complete.obs")*cor(y,z, use="complete.obs"))/
                                     sqrt((1-cor(x,z, use="complete.obs")^2)*(1-cor(y,z, use="complete.obs")^2)),4)}

##(1) ANALYSIS FOR FOUR CORPORA
 
## (1) MRC (Wilson, 1988)
norms = read.csv("MRC_corpus.csv",header=TRUE)
freqs = read.table("SUBTLEXusDataBase.txt",header=TRUE)

# --add frequency data --
index <- match(norms$word, freqs$Word)
norms$freq <- freqs$Lg10WF[index]

#concreteness
c_E_cp = cor.test(norms$mrc.phon,norms$mrc.conc, use="complete.obs")
pc_E_cp = bm.partial(norms$mrc.phon,norms$mrc.conc,norms$freq)

c_E_cl = cor.test(norms$mrc.nlet,norms$mrc.conc, use="complete.obs")
pc_E_cl = bm.partial(norms$mrc.nlet,norms$mrc.conc,norms$freq)

c_E_cs = cor.test(norms$mrc.syl,norms$mrc.conc, use="complete.obs")
pc_E_cs = bm.partial(norms$mrc.syl,norms$mrc.conc,norms$freq)

summary(lm(norms$mrc.nlet ~ norms$mrc.conc + norms$freq))

#familiarity
c_E_fp = cor.test(norms$mrc.phon,norms$mrc.fam, use="complete.obs")
pc_E_fp = bm.partial(norms$mrc.phon,norms$mrc.fam,norms$freq)

c_E_fl = cor.test(norms$mrc.nlet,norms$mrc.fam, use="complete.obs")
pc_E_fl = bm.partial(norms$mrc.nlet,norms$mrc.fam,norms$freq)

c_E_fs = cor.test(norms$mrc.syl,norms$mrc.fam, use="complete.obs")
pc_E_fs =bm.partial(norms$mrc.syl,norms$mrc.fam,norms$freq)

summary(lm(norms$mrc.nlet ~ norms$mrc.fam))
summary(lm(norms$mrc.nlet ~ norms$mrc.fam + norms$freq))

#imageability
c_E_ip = cor.test(norms$mrc.phon,norms$mrc.imag, use="complete.obs")
pc_E_ip = bm.partial(norms$mrc.phon,norms$mrc.imag,norms$freq)

c_E_il = cor.test(norms$mrc.nlet,norms$mrc.imag, use="complete.obs")
pc_E_il = bm.partial(norms$mrc.nlet,norms$mrc.imag,norms$freq)

c_E_is = cor.test(norms$mrc.syl,norms$mrc.imag, use="complete.obs")
pc_E_is = bm.partial(norms$mrc.syl,norms$mrc.imag,norms$freq)

summary(lm(norms$mrc.nlet ~ norms$mrc.imag + norms$freq))

# correlation between frequency and familiartiy
cor.test(norms$mrc.fam,norms$freq, use="complete.obs")
cor.test(norms$mrc.conc,norms$freq, use="complete.obs")
cor.test(norms$mrc.imag,norms$freq, use="complete.obs")


## (2) Italian (Della Rosa et al., 2010) 
di = read.csv("DellaRosaEtAl_norms.csv",header=TRUE)
di = di[,1:17] #remove extra colums

freqsI =  read.csv("lemmi-export.csv", header=TRUE)
freqsC = aggregate(freqsI$linked.tot, by=list(freqsI$lemma), FUN=sum) # collapse across different word classes
names(freqsC) <- c("word", "freq")
freqsC$log_freq = log(freqsC$freq)

# --add frequency stats --
di$Italian.name <- tolower(di$Italian.name)
index <- match(di$Italian.name, freqsC$word)
di$freq <- freqsC$log_freq[index]

#concretness
c_I_cl = cor.test(di$CNC,di$Let)
pc_I_cl = bm.partial(di$Let,di$CNC,di$freq)

#imageability
c_I_il = cor.test(di$IMG,di$Let)
pc_I_il = bm.partial(di$Let,di$IMG,di$freq)

#familiarty
c_I_fl = cor.test(di$FAM,di$Let)
pc_I_fl = bm.partial(di$Let,di$FAM,di$freq)

#abstractness
c_I_al = cor.test(di$ABS,di$Let)
pc_I_al = bm.partial(di$Let,di$ABS,di$freq)

summary(lm(di$Let ~ di$CNC + di$freq))
summary(lm(di$Let ~ di$IMG + di$freq))
summary(lm(di$Let ~ di$FAM + di$freq))
summary(lm(di$Let ~ di$ABS + di$freq))


## (3) French (Desrochers & Thompson, 2009)
d = read.csv("Desrochers-Thompson_2009_Ratings.csv",header=TRUE, stringsAsFactors=FALSE, encoding = "UTF-8")

# get length
d$length = nchar(d$NOUN)

c_F_il = cor.test(d$IMAGE_Mean,d$length)
pc_F_il = bm.partial(d$length, d$IMAGE_Mean, log(d$FREQ_Mean))

summary(lm(d$length ~ d$IMAGE_Mean+ log(d$FREQ_Mean)))


## (4) Concretness  (Brysbaert et al., 2013)
b <- read.csv("brysbaert_corpus.csv",header=TRUE)
b <- b[b$Word != "",] # get rid of empty rows
b <- b[b$Bigram == 0,]# get rid of two word lemmas

# --add frequency stats --
index <- match(b$Word, freqs$Word)
b$logfreq <- freqs$Lg10WF[index]

c_C_cl = cor.test(b$Length,b$Conc.M)
pc_C_cl = bm.partial(b$Length,b$Conc.M,b$logfreq)

summary(lm(b$Length ~ b$Conc.M))
summary(lm(b$Length ~ b$Conc.M + b$logfreq))



##(2) MAKE PLOT
colors = brewer.pal(7, "Set1")
use_colors = c(colors[1], colors[2], colors[3], colors[4])
conc_color = use_colors[1]
fam_color = use_colors[2]
img_color = use_colors[3]
abs_color = use_colors[4]
small_space = .2
big_space = 1.4
corpus_position = -.48
corpus_l_postion = -.51
lwd_weight = 3
os = .4

correlations =  c( c_E_cl$estimate,     c_E_fl$estimate,    c_E_il$estimate, c_C_cl$estimate, c_I_cl$estimate, c_I_fl$estimate, c_I_il$estimate,    c_I_al$estimate, c_F_il$estimate)  
cihigh =   c( c_E_cl$conf.int[1],  c_E_fl$conf.int[1],  c_E_il$conf.int[1], c_C_cl$conf.int[1], c_I_cl$conf.int[1], c_I_fl$conf.int[1], c_I_il$conf.int[1],    c_I_al$conf.int[1], c_F_il$conf.int[1])  
cilow =   c( c_E_cl$conf.int[2],  c_E_fl$conf.int[2],  c_E_il$conf.int[2], c_C_cl$conf.int[2], c_I_cl$conf.int[2], c_I_fl$conf.int[2], c_I_il$conf.int[2],    c_I_al$conf.int[2], c_F_il$conf.int[2])  
partial_correlations = c( pc_E_cl,    pc_E_fl,   pc_E_il, pc_I_cl,pc_C_cl, pc_I_fl,   pc_I_il,   pc_I_al, pc_F_il)
measure_colors =       c( conc_color, fam_color, img_color, conc_color, conc_color, fam_color, img_color, abs_color,  img_color)

bplt = barplot(correlations, 
        ylab="Correlation with length",
        col=measure_colors,
        ylim = c(-.5, .5), 
        names.arg = c(""),
        space = c(small_space, small_space, small_space, big_space, big_space, small_space, 
                         small_space, small_space, big_space)
)
arrows(bplt,cilow, bplt,cihigh, angle=90, code=3, length=0) #cis
text(x= bplt, y=  partial_correlations, labels=as.character('.'), xpd=TRUE, cex = 5)
axis(2, at=seq(-.5,.5,.1))
abline(0, 0)
legend(.2, .48, c("concreteness", 
                  "familiarity", 
                  "imageability", 
                  "abstractness"), 
       fill= use_colors, bty = "n")

# add corpus labels
segments(bplt[1]-os, corpus_position, bplt[3]+os, corpus_position, lwd = lwd_weight )
segments(bplt[4]-os, corpus_position, bplt[4]+os, corpus_position, lwd = lwd_weight )
segments(bplt[5]-os, corpus_position, bplt[8]+os, corpus_position, lwd = lwd_weight )
segments(bplt[9]-os, corpus_position, bplt[9]+os, corpus_position, lwd = lwd_weight )

corpus_l_postion = -.51
text(x=bplt[2],corpus_l_postion, labels=as.character('English'), xpd=TRUE)
text(x=bplt[2],corpus_l_postion -.03, labels=as.character('(MRC)'), xpd=TRUE)
text(x=bplt[4],corpus_l_postion, labels=as.character('English'), xpd=TRUE)
text(x=bplt[4],corpus_l_postion -.03, labels=as.character('(Brysbaert et al., 2013)'), xpd=TRUE)
text(x=(bplt[6] + bplt[7])/2,corpus_l_postion, labels=as.character('Italian'), xpd=TRUE)
text(x=bplt[9],corpus_l_postion, labels=as.character('French'), xpd=TRUE)
