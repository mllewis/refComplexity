# Experimental data analysis for Cogsci2014 Referential Complexity paper (Lewis & Frank)
#
# Data for 3 adult experiments are saved as 1 file (refComplex_adult.csv), 
# and data for kid experiment as a second file (refComplex_kid.csv).
# Experiment labels do not correspond to labels in the paper. Here is the key:
# Exp. 1a=4 (geons, fc), Exp. 1b=3 (geons, bet), Exp. 2=12 (real objects), 15 (kids #3)
#
# This script draws the plots presented in the paper (part 1), and reproduces the analysis (part 2).

rm(list=ls())

# load packages
library(reshape)
library(plyr)
library(bootES)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(bootstrap)
library(lme4)

# define functions
getSummary <- function(D){
  summary_fc1 <- ddply(D[D$measure == "fc",], .(exp, condition), function (d) {   
    b <- boot(d$response, function(u,i) table(u[i])["complex"]/length(u), R = 1000) 
    ci <- boot.ci(b, type =  "basic")
    ciwl = ci[4][[1]][1,4]
    ciwu = ci[4][[1]][1,5]
    return (c(ciwl, ciwu))})
  names(summary_fc1)[3:4] = c("cill", "ciul")            
  summary_fc2 <- ddply(D[D$measure == "fc",], .(exp, condition), summarize, 
                       complex_proport = sum(response=="complex")/length(response),
                       n = length(workerid), 
                       exp_date = creationtime[1])
  summary_fc <- merge(summary_fc1, summary_fc2)             
  summary_bet1 <- ddply(D[D$measure == "bet",], .(exp, condition), function(d) {(sd(as.numeric(as.character(d$response))) / sqrt(length(as.numeric(as.character(d$response))))) * .0196})
  names(summary_bet1)[3]= "err"            
  summary_bet2 <- ddply(D[D$measure == "bet",], .(exp, condition), summarize, 
                        complex_proport = mean(as.numeric(as.character(response)))/100, 
                        n = length(workerid), 
                        exp_date = creationtime[1],
                        sd = sd(response)/100)
  summary_bet <- merge(summary_bet1, summary_bet2) 
  summary_bet$cill <- summary_bet$complex_proport - summary_bet$err   
  summary_bet$ciul <- summary_bet$complex_proport + summary_bet$err
  #summary_bet$err <- NULL                          
  metaD <- rbind.fill(summary_bet, summary_fc)
  metaD <- metaD[!is.na(metaD$condition),] # remove empty exps
  return(metaD)
} 
getSummary_other <- function(D){
  summary_fc1 <- ddply(D[D$measure == "fc",], .(exp, condition_other), function (d) {   
    b <- boot(d$response, function(u,i) table(u[i])["complex"]/length(u), R = 1000) 
    ci <- boot.ci(b, type =  "basic")
    ciwl = ci[4][[1]][1,4]
    ciwu = ci[4][[1]][1,5]
    return (c(ciwl, ciwu))})
  names(summary_fc1)[3:4] = c("cill", "ciul")            
  summary_fc2 <- ddply(D[D$measure == "fc",], .(exp, condition_other), summarize, 
                       complex_proport = sum(response=="complex")/length(response),
                       n = length(workerid), 
                       exp_date = creationtime[1])
  summary_fc <- merge(summary_fc1, summary_fc2)             
  summary_bet1 <- ddply(D[D$measure == "bet",], .(exp, condition_other), function(d) {(sd(as.numeric(as.character(d$response))) / sqrt(length(as.numeric(as.character(d$response))))) * .0196})
  names(summary_bet1)[3]= "err"            
  summary_bet2 <- ddply(D[D$measure == "bet",], .(exp, condition_other), summarize, 
                        complex_proport = mean(as.numeric(as.character(response)))/100, 
                        n = length(workerid), 
                        exp_date = creationtime[1],
                        sd = sd(response)/100)
  summary_bet <- merge(summary_bet1, summary_bet2) 
  summary_bet$cill <- summary_bet$complex_proport - summary_bet$err   
  summary_bet$ciul <- summary_bet$complex_proport + summary_bet$err
  summary_bet$err <- NULL                          
  metaD <- rbind.fill(summary_bet, summary_fc)
  metaD <- metaD[!is.na(metaD$condition_other),] # remove empty exps
  return(metaD)
} 
d.fc <- function(d) {
  cond <- all(intersect(levels(d$condition),c("long","short")) == c("long","short"))
  
  if (cond) {
    d<- d[d$condition == "long" | d$condition == "short",]
    d <- droplevels(d)
    ns = table(d$condition, d$response)
    or = (ns[1]*ns[4])/(ns[2]*ns[3])
    cf = pi/sqrt(3)
    effect_size = log(or)/cf
    se = sqrt(1/(ns[1] + 1/ns[2] + 1/ns[3] + 1/ns[4]))
    d_se = se * (3/(pi^2))
    d_err = d_se*1.96
    cill = effect_size - d_err
    ciul = effect_size + d_err
    es <- data.frame(effect_size=effect_size,
                     cill = cill,
                     ciul = ciul)
  } else {
    es <- data.frame(effect_size=NA,
                     cill = NA,
                     ciul = NA)
  }
  
  return (es)
}
d.bet <- function(d){
  d$response.num <- as.numeric(as.character(d$response))
  
  cond <- mean(d$response.num)>1 &&
    all(sort(drop.levels(d$condition)[1:2]) == c("long" , "short"), na.rm = TRUE) &&
    nrow(d) > 1
  
  if (!is.na(cond) & cond) {
    boots <- bootES(d, 2000,
                    data.col = "response.num",
                    group.col = "condition",
                    contrast = c("short", "long"),
                    effect.type =  "cohens.d")
    
    es <- data.frame(effect_size=boots$t0,
                     cill = boots$bounds[1],
                     ciul = boots$bounds[2])
    
  } else {
    es <- data.frame(effect_size=NA,
                     cill = NA,
                     ciul = NA)
    
  }
  
  return(es)
}
theta <- function(x,xdata,na.rm=T) {mean(xdata[x],na.rm=na.rm)}
ci.low <- function(x,na.rm=T) {
  mean(x,na.rm=na.rm) - quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.025,na.rm=na.rm)}
ci.high <- function(x,na.rm=T) {
  quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.975,na.rm=na.rm) - mean(x,na.rm=na.rm)}

# load data
D <- read.csv("refComplex_adult.csv")
D_good <- D[D$duplicateS== "no" & (D$filter == TRUE | is.na(D$filter)),] # exclude bad check and duplicate subjects
d <- read.csv("refComplex_kid.csv")

### (1) DRAW PLOTS
# parameters
small_space = .2
big_space = .4
huge_space = 1.1
h1 = -.03
h2 = -.09
h3 = -.15
h4 = -.19
lwd_weight = 3
os = .3

# set up colors
colors = brewer.pal(7, "Set1")
use_colors = c(colors[1], colors[2], colors[3], colors[4])
measure_colors = c(use_colors[2], use_colors[1])
control_blue = "#c3d8e9"
control_red = "#f6baba"

# PREP DATA
# E1a
E1as= getSummary(D_good[D_good$exp2 == 4,])[2,"complex_proport"]
E1al= getSummary(D_good[D_good$exp2 == 4,])[1,"complex_proport"]
E1as_l= getSummary(D_good[D_good$exp2 == 4,])[2,"cill"]
E1al_l= getSummary(D_good[D_good$exp2 == 4,])[1,"cill"]
E1as_h= getSummary(D_good[D_good$exp2 == 4,])[2,"ciul"]
E1al_h= getSummary(D_good[D_good$exp2 == 4,])[1,"ciul"]

# E1b 
d3 <-  getSummary_other(D_good[D_good$exp2 == 3,])
d3$language <- c("long", "long", "long", "short", "short",  "short")
d3$trial.type <- c("complex-complex", "simple-complex", "simple-simple", 
                   "complex-complex", "simple-complex", "simple-simple")
E1bscc= d3[d3$trial.type == "complex-complex" & d3$language == "short","complex_proport"]
E1blcc= d3[d3$trial.type == "complex-complex" & d3$language == "long","complex_proport"]
E1bssc= d3[d3$trial.type == "simple-complex" & d3$language == "short","complex_proport"]
E1blsc= d3[d3$trial.type == "simple-complex" & d3$language == "long","complex_proport"]
E1bsss= d3[d3$trial.type == "simple-simple" & d3$language == "short","complex_proport"]
E1blss= d3[d3$trial.type == "simple-simple" & d3$language == "long","complex_proport"]
E1bscc_l = d3[d3$trial.type == "complex-complex" & d3$language == "short","cill"]
E1blcc_l = d3[d3$trial.type == "complex-complex" & d3$language == "long","cill"]
E1bssc_l = d3[d3$trial.type == "simple-complex" & d3$language == "short","cill"]
E1blsc_l = d3[d3$trial.type == "simple-complex" & d3$language == "long","cill"]
E1bsss_l = d3[d3$trial.type == "simple-simple" & d3$language == "short","cill"]
E1blss_l = d3[d3$trial.type == "simple-simple" & d3$language == "long","cill"]
E1bscc_h = d3[d3$trial.type == "complex-complex" & d3$language == "short","ciul"]
E1blcc_h = d3[d3$trial.type == "complex-complex" & d3$language == "long","ciul"]
E1bssc_h = d3[d3$trial.type == "simple-complex" & d3$language == "short","ciul"]
E1blsc_h = d3[d3$trial.type == "simple-complex" & d3$language == "long","ciul"]
E1bsss_h = d3[d3$trial.type == "simple-simple" & d3$language == "short","ciul"]
E1blss_h = d3[d3$trial.type == "simple-simple" & d3$language == "long","ciul"]

# E2
E2s=getSummary(D_good[D_good$exp2 == 12,])[2,"complex_proport"]
E2l=getSummary(D_good[D_good$exp2 == 12,])[1,"complex_proport"]
E2s_l= getSummary(D_good[D_good$exp2 == 12,])[2,"cill"]
E2l_l= getSummary(D_good[D_good$exp2 == 12,])[1,"cill"]
E2s_h=getSummary(D_good[D_good$exp2 == 12,])[2,"ciul"]
E2l_h=getSummary(D_good[D_good$exp2 == 12,])[1,"ciul"]

#E3
d <- d[which(d$trunc_age != 2),] #throw out two year olds
d <- d[which(d$exclusion == "-"),] #throw out trials to be excluded
d$male <- as.factor(d$male)
sub.d <- d[d$critORfill == "critical",] #critical only 
sub.d$propComplex <- sub.d$response == "C" #make response boolean

#compute mean correct by kid
mss <- aggregate(propComplex ~ trialtype + trunc_age + subid, data=sub.d,FUN=mean)

#compute mean correct across kids
d15 <- aggregate(propComplex ~ trialtype + trunc_age, data=mss,FUN=mean)
d15$cih <- aggregate(propComplex ~ trialtype + trunc_age, data=mss,FUN=ci.high)$propComplex
d15$cil <- aggregate(propComplex ~ trialtype + trunc_age, data=mss,FUN=ci.low)$propComplex
levels(d15$trialtype) <- c("long", "short")

E4_3_s=d15[2, "propComplex"]
E4_3_l=d15[1, "propComplex"]
E4_4_s =d15[4, "propComplex"]
E4_4_l=d15[3, "propComplex"]
E4_5_s=d15[6, "propComplex"]
E4_5_l=d15[5, "propComplex"]
E4_3_s_l= E4_3_s - d15[2, "cil"] 
E4_3_l_l= E4_3_l - d15[1, "cil"] 
E4_4_s_l= E4_4_s - d15[4, "cil"] 
E4_4_l_l= E4_4_l - d15[3, "cil"] 
E4_5_s_l= E4_5_s - d15[6, "cil"] 
E4_5_l_l= E4_5_l - d15[5, "cil"] 
E4_3_s_h=d15[2, "cih"] + E4_3_s
E4_3_l_h=d15[1, "cih"] + E4_3_l
E4_4_s_h=d15[4, "cih"] + E4_4_s
E4_4_l_h=d15[3, "cih"] + E4_4_l
E4_5_s_h=d15[6, "cih"] + E4_5_s
E4_5_l_h=d15[5, "cih"] + E4_5_l

#merge experimental data together
exp_meansA =  c( E1as, E1al, E1bscc, E1blcc, E1bssc, E1blsc,E1bsss, E1blss, E2s, E2l)
exp_low_ciA =  c(E1as_l,E1al_l, E1bscc_l, E1blcc_l, E1bssc_l,E1blsc_l,E1bsss_l, E1blss_l, E2s_l, E2l_l)
exp_hi_ciA =  c(E1as_h,E1al_h,  E1bscc_h, E1blcc_h, E1bssc_h,E1blsc_h,E1bsss_h, E1blss_h, E2s_h, E2l_h)
bar_colorsA =  c( measure_colors[1], measure_colors[2], control_blue, control_red,measure_colors[1], measure_colors[2], control_blue, control_red,measure_colors[1], measure_colors[2])
exp_meansB =  c( E4_3_s, E4_3_l, E4_4_s, E4_4_l, E4_5_s, E4_5_l)
exp_low_ciB =  c( E4_3_s_l, E4_3_l_l, E4_4_s_l, E4_4_l_l, E4_5_s_l, E4_5_l_l)
exp_hi_ciB =  c(E4_3_s_h, E4_3_l_h, E4_4_s_h, E4_4_l_h, E4_5_s_h, E4_5_l_h)
bar_colorsB =  c(measure_colors[1], measure_colors[2],measure_colors[1], measure_colors[2],measure_colors[1], measure_colors[2])

#Draw plot
#plot A
brpltA = barplot(exp_meansA, 
               ylab="Proportion generalizations to complex object",
               col=bar_colorsA,
               ylim = c(0,1), 
               space = c(small_space, small_space, huge_space, small_space, big_space,small_space, 
                         big_space, small_space, huge_space, small_space),
               axes= TRUE, 
               cex.lab=1.3
)
arrows(brpltA,exp_low_ciA, brpltA,exp_hi_ciA, angle=90, code=3, length=0) #cis
abline(a=.5,b= 0, col="black", lty=2) #chance line

# experiment headings
text((brpltA[3] + brpltA[4])/2,h1, labels=as.character('comp./comp.'), xpd=TRUE, cex = 1.1)
text((brpltA[5] + brpltA[6])/2,h1, labels=as.character('comp./simp.'), xpd=TRUE, cex = 1.1)
text((brpltA[7] + brpltA[8])/2,h1, labels=as.character('simp./simp.'), xpd=TRUE, cex = 1.1)
text((brpltA[3] + brpltA[8])/2,h2, labels=as.character('Objects'), xpd=TRUE, cex = 1.2)

# experiment num labels
segments(brpltA[1]-os, h3, brpltA[2] +os, h3, lwd = lwd_weight, xpd = TRUE)
segments(brpltA[3]-os, h3, brpltA[8] +os, h3, lwd = lwd_weight, xpd = TRUE)
segments(brpltA[9]-os, h3, brpltA[10] +os, h3, lwd = lwd_weight, xpd = TRUE)
text((brpltA[1] + brpltA[2])/2,h4, labels=as.character('Exp. 1a'), xpd=TRUE, cex = 1.3)
text((brpltA[3] + brpltA[8])/2,h4, labels=as.character('Exp. 1b'), xpd=TRUE, cex = 1.3)
text((brpltA[9] + brpltA[10])/2,h4, labels=as.character('Exp. 2'), xpd=TRUE, cex = 1.3)

#plot B
brpltB = barplot(exp_meansB, 
                 ylab="Proportion generalizations to complex object",
                 col=bar_colorsB,
                 ylim = c(0,1), 
                 space = c(huge_space,  small_space, big_space,  small_space, big_space, small_space),
                 axes= TRUE, 
                 cex.lab=1.3
)
arrows(brpltB,exp_low_ciB, brpltB,exp_hi_ciB, angle=90, code=3, length=0) #cis
abline(a=.5,b= 0, col="black", lty=2) #chance line
text((brpltB[1] + brpltB[2])/2,h1, labels=as.character('3'), xpd=TRUE, cex = 1.1)
text((brpltB[3] + brpltB[4])/2,h1, labels=as.character('4'), xpd=TRUE, cex = 1.1)
text((brpltB[5] + brpltB[6])/2,h1, labels=as.character('5'), xpd=TRUE, cex = 1.1)
text((brpltB[1] + brpltB[6])/2,h2, labels=as.character('Age (years)'), xpd=TRUE, cex = 1.2)

segments(brpltB[1]-os, h3, brpltB[6] +os, h3, lwd = lwd_weight, xpd = TRUE)
text((brpltB[1] + brpltB[6])/2,h4, labels=as.character('Exp. 3'), xpd=TRUE, cex = 1.3)


### (2) STATS
#Exp 4 (1a)
d4l = D_good[D_good$exp2 == 4,]
d4l$response = drop.levels(d4l$response)
d4l$condition = drop.levels(d4l$condition)

## compare long and short
chisq.test(table(d4l$response, d4l$condition))

#<<compare without exclusions>>
d4lBAD = D[D$exp2 == 4,]
d4lBAD$response = drop.levels(d4lBAD$response)
d4lBAD$condition = drop.levels(d4lBAD$condition)

## compare long and short
chisq.test(table(d4lBAD$response, d4lBAD$condition))

## get d
d.fc(d4l)

#Exp 3 (1b)
d3l = D_good[D_good$exp2 == 3,]

#add object alternatives variable
d3l$object <- ifelse(d3l$condition_other == '\"complex_sc\"' | d3l$condition_other == '\"simple_sc\"', "sc", 
                     ifelse(d3l$condition_other == '\"complex_cc\"' | d3l$condition_other == '\"simple_cc\"', "cc",
                            'ss'))
d3l$object  = as.factor(d3l$object)

#compare long vs. short for each alternatives condition
t.test(as.numeric(as.character(response))~condition,d3l[d3l$object == "sc",], var.equal= TRUE)
t.test(as.numeric(as.character(response))~condition_other,d3l[d3l$object == "ss",], var.equal= TRUE)
t.test(as.numeric(as.character(response))~condition_other,d3l[d3l$object == "cc",], var.equal= TRUE)

#<<compare without exclusions>>
d3lBAD = D[D$exp2 == 3,]
#add object alternatives variable
d3lBAD$object <- ifelse(d3lBAD$condition_other == '\"complex_sc\"' | d3lBAD$condition_other == '\"simple_sc\"', "sc", 
                     ifelse(d3lBAD$condition_other == '\"complex_cc\"' | d3lBAD$condition_other == '\"simple_cc\"', "cc",
                            'ss'))
d3lBAD$object  = as.factor(d3lBAD$object )

#compare long vs. short for each alternatives condition
t.test(as.numeric(as.character(response))~condition,d3lBAD[d3lBAD$object == "sc",], var.equal= TRUE)
t.test(as.numeric(as.character(response))~condition_other,d3lBAD[d3lBAD$object == "ss",], var.equal= TRUE)
t.test(as.numeric(as.character(response))~condition_other,d3lBAD[d3lBAD$object == "cc",], var.equal= TRUE)

#get effect size
d.bet(d3l[d3l$object == "sc",])

#Exp 12 (2)
# get exclusions counts
d12l = D_good[D_good$exp2 == 12,]

#compare long to short
t.test(as.numeric(as.character(response))~condition,d12l, var.equal= TRUE)

#<<compare without exclusions>>
d12lBAD = D[D$exp2 == 12,]
t.test(as.numeric(as.character(response))~condition,d12lBAD, var.equal= TRUE)

#get effect size
d.bet(d12l)

# Exp 15 (3)
lmer(response ~ trialtype + (1|subid) + (trunc_age|trialnum),
     family="binomial",sub.d )

lmer(response ~ trunc_age*trialtype + (1|subid) + (trunc_age|trialnum),
     family="binomial",sub.d  )

lmer(response ~ trunc_age*trialtype + (1|subid) + (1|trialnum),
     family="binomial",sub.d )

lmer(response ~ age*trialtype + (1|subid) + (age|trialnum), #The model reported in the paper
     family="binomial",sub.d  )

lmer(response ~ age*trialtype + (1|subid) + (1|trialnum),
     family="binomial",sub.d  )

# age groups
t.test(propComplex~trialtype,mss[mss$trunc_age==3,], paired = TRUE)
t.test(propComplex~trialtype,mss[mss$trunc_age==4,], paired = TRUE)
t.test(propComplex~trialtype,mss[mss$trunc_age==5,], paired = TRUE)

#get effect size
bootES(mss[mss$trunc_age == 5,], 2000,
                data.col = "propComplex",
                group.col = "trialtype",
                contrast = c("S", "L"),
                effect.type =  "cohens.d")
