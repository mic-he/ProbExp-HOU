### preliminary data cleaning and exploration for the complex expressions experiments

## libraries, sources, data
# manipulation
library(dplyr)
library(reshape2)
# visualization
library(ggplot2)
library(RColorBrewer)
library(hexbin)
library(gridExtra)
# analysis
library(mlogit) # multinomial logistic regression
library(VGAM) # dbetabinom.ab distribution is here
library(bootstrap) # boostrap 95% CIs

# source functions
source("useful_functions.R") # various useful definitions
source("prediction_functions.R") # the prediction functions, ie R implementation of the model 

# get data
df.production <- read.csv("raw_production_data.csv") # raw data from complex expression production task
df.interpretation <- read.csv("raw_interpretation_data.csv") # raw data from complex expression interpretation task
df.observed.belief <- read.csv("observed_belief.csv") # processed data from belief-about-the-urn task

## production task: clean and explore
xdf.production <- df.production

# how many subjects?
xdf.production$id %>% unique() %>% length()

# randomization of experimental conditions
xtabs(~value, xdf.production)

# reported languages
xdf.production$language %>% levels()
# [1] "American English" "egnlish"          "Engligh"          "english"          "eNGLISH"          "English"          "ENGLISH"
# spelling issues aside, no problems here

# comments
xdf.production$comments %>% levels()
# interesting output:
# [3] "Always a 50/50 chance unless all the balls are drawn and are all the same color"
# [7] "I accidentally responded the opposite to one question because I was thinking about the blue rather than the red."
# [10] "I would have liked to have a couple more options in terms of possible answers, but overall a very clear and easy task."
# [27] "The choices for messages to send were rather vague and didn't translate numbers well"
# [28] "The responses felt very structured.  In the real world I probably wouldn't phrase things like that. "

# we might want to exclude data from subjects corresponding to comments 3 and 7
ids=droplevels(subset(xdf.production, xdf.production$comments=="Always a 50/50 chance unless all the balls are drawn and are all the same color" 
                      | xdf.production$comments=="I accidentally responded the opposite to one question because I was thinking about the blue rather than the red."))$id
xdf.production=droplevels(subset(xdf.production, !xdf.production$id %in% ids ))

# how many subjects?
xdf.production$id %>% unique() %>% length()
# 102

# reported difficulty
xtabs(~difficulty, xdf.production)/12 # 12 is the number of trials completed by each participant
# difficulty
# 0   1   2   3   4   5   6   7   8  10 # levels of difficulty  
# 15  12  28  13  11  11  2   5   4 1   # counts
# no problems here

# group inner occurrences of likely and probable
xdf.production$answer3=ifelse(xdf.production$answer2=="probable" | xdf.production$answer2=="likely","likely",levels(xdf.production$answer2)[xdf.production$answer2])

# produce and save a more compact, better looking data frame with less "useless" columns and a couple more useful columns
xdf.production$comments<-xdf.production$engagement<-xdf.production$difficulty<-xdf.production$language<-NULL
# add column with full complex answer
xdf.production$answer.full=factor(paste(xdf.production$answer1,xdf.production$answer3))
# add column with level of high-order uncertainty
xdf.production$uncertainty=ifelse(xdf.production$access==10,"none",ifelse(xdf.production$access>5,"low","high"))
# save (to be used in analysis.r)
write.csv(xdf.production,"clean_production_data.csv", row.names = FALSE)

# visualize expression choices in each condition with bar plots
# get counts of expression per condition
counts=t(xtabs(~answer.full+value,xdf.production))
# how many observation for each condition?
observations=as.vector(xtabs(~value,xdf.production))
# counts to percentage, as data frame
df.e=as.data.frame(melt(counts/observations*100))
colnames(df.e)=c("value","expression","percentage")
# add some columns to df for easier visualization
df.e$access=trunc(df.e$value/10)
df.e$observation=df.e$value-10*df.e$access
# better looking labels for values
df.e$label=as.factor(paste0(df.e$observation," red balls out of ",df.e$access))
# change order to make plot more readable
df.e$label = factor(df.e$label,levels(df.e$label)[c(1,5,3,6,9,2,7,10,12,15,4,8,11,13,14)])
# better looking names for expressions
levels(df.e$expression)=c("is likely","is possible", "is unlikely", "is certainly likely", "is certainly possible", "is certainly unlikely", "is probably likely", "is probably possible", "is probably unlikely", "might be likely", "might be possible", "might be unlikely")
# change order to make plot more readable
df.e$expression = factor(df.e$expression,levels(df.e$expression)[c(12,11,10,9,8,7,6,5,4,3,2,1)])
# plot
expression.bars=ggplot(data=df.e)+geom_bar(aes(x=expression,y=percentage,fill=expression),stat="identity")+ scale_y_continuous(limits = c(0,60), breaks=seq(0,60,10))+coord_flip()+facet_wrap(facets = ~label, ncol=5)+theme_bw()+theme(strip.background = element_blank(),text=element_text(size=24), legend.position="none")
show(expression.bars)


## add bootstrapped 95% CIs
values <- levels(as.factor(xdf.production$value)) # experimental conditions
expressions <- levels(xdf.production$answer.full)
repetitions <- 1000 # how many samples?

# output is array values x expressions x repetitions, output1 is values x expressions x 3 (ci.low, mean, ci.high)
output <- array(dim=c(length(values),length(expressions),repetitions),dimnames = list(values,expressions,c(1:repetitions)))
output1 <- array(dim=c(length(values),length(expressions),3),dimnames = list(values,expressions,c("ci.low","mean","ci.high")))
for (v in 1:length(values)){ # for each condition...
  choices <- subset(xdf.production, xdf.production$value==values[v])$answer.full # get the vector with observed answers
  for (e in 1:length(expressions)){
    output[v,e,] <- bootstrap(x=choices,nboot=repetitions,theta = percentage,expressions[e])$thetastar
    output1[v,e,] <- c(ci.low(output[v,e,]),mean(output[v,e,]),ci.high(output[v,e,]))
  }
}
# as df
df.e.cis <- as.data.frame(cbind(melt(output1[,,1]),melt(output1[,,2])$value,melt(output1[,,3])$value))
colnames(df.e.cis) <- c("value","message","ci.low","mean","ci.high")
# save
write.csv(df.e.cis,"complex_production_cis.csv", row.names = FALSE)



############################################################################


## interpretation task: clean and explore
xdf.interpretation <- df.interpretation

# how many subjects?
xdf.interpretation$id %>% unique() %>% length()
#158

# reported languages
xdf.interpretation$language %>% levels()
#[1] "engish"     "Engish"     "english"    "English"    "ENGLISH"    "english "   "hindi"      "Macedonian" "Spanish"    "Tamil"     
#[11] "Telugu" 

# drop non native english speakers
xdf.interpretation=droplevels(subset(xdf.interpretation,xdf.interpretation$language!="Macedonian" &
                                       xdf.interpretation$language!="Tamil" &
                                       xdf.interpretation$language!="Telugu" &
                                       xdf.interpretation$language!="hindi" &
                                       xdf.interpretation$language!="Spanish"))
# how many subjects left?
xdf.interpretation$id %>% unique() %>% length()
#151

# comments
xdf.interpretation$comments %>% levels()
# interesting output:
#[5] "I don't know what you're studying, but I acted upon the suggestions of \"possible\" and \"not possible\" and \"maybe\" and etc. I hope I did this right. I hope there is not a correct or an incorrect answer, in which case I would have to just give up in confusion."
#[6] "I don't think the experiment made sense to me and I'm not 100% sure of what I just did. How would a person be able to guess how many red balls are in a jar based on a vague statement? No one would be accurate unless they were very lucky."                          
#[7] "I felt like you can't make any assumptions about the balls that you haven't the odds for each ball are the same regardless of balls already pulled."                                                                                                                    
#[8] "I know the task was worded this way specifically for the scenario, but it's unlikely that anyone would actually speak in this manner in reality, as double negatives have a tendency to confuse listeners."    
#[10] "i think I screwed up the math quite a bit"                                                                                                                                                                                                                              
#[11] "I think it would be nice to see both sides of the interaction during the practice round so that you understand it more fully." 
#[14] "It was easy-ish to guess but very difficult to be right.  It was fun though, thanks!"  
#[17] "It was interesting. It reminded me how sords can be tricky sometimes."                                                                                                                                                                                                  
#[18] "It was somewhat difficult to reason out the statements"  
#[29] "No, it was quite difficult however it was an overall good experience because you need to think from your partner's perspective." 

# reported difficulty
xtabs(~difficulty, xdf.interpretation)/24 # 24 is the number of trials completed by each participant
# difficulty
#0  1  2  3  4  5  6  7  8  9 10 
#26 14 24 16 14 23 16 11  8  2  2
#ok

# has anybody said acc=0 and obs=0 for any expression?
xdf.interpretation %>% filter(access==0 & observation==0) %>% select(id) %>% n_distinct()
# 1
# drop it
pessimists=levels(droplevels(subset(xdf.interpretation, xdf.interpretation$access==0 & xdf.interpretation$observation==0))$id)
xdf.interpretation=droplevels(subset(xdf.interpretation,!xdf.interpretation$id %in% pessimists))
# how many subjects left?
xdf.interpretation$id %>% unique() %>% length()
#150

# produce and save a more compact, better looking data frame
# rename expression levels, to make them more readable and match names in production experiment and model
xdf.interpretation$expression=ifelse(xdf.interpretation$expression=="cerLik","is_certainly likely",
                              ifelse(xdf.interpretation$expression=="cerPos","is_certainly possible",
                              ifelse(xdf.interpretation$expression=="cerUnlik","is_certainly unlikely",
                              ifelse(xdf.interpretation$expression=="isLik","is likely",
                              ifelse(xdf.interpretation$expression=="isPos","is possible",
                              ifelse(xdf.interpretation$expression=="isUnlik","is unlikely",
                              ifelse(xdf.interpretation$expression=="mightLik","might_be likely",
                              ifelse(xdf.interpretation$expression=="mightPos","might_be possible",
                              ifelse(xdf.interpretation$expression=="mightUnlik","might_be unlikely",
                              ifelse(xdf.interpretation$expression=="probLik","is_probably likely",
                              ifelse(xdf.interpretation$expression=="probPos","is_probably possible",
                              ifelse(xdf.interpretation$expression=="probUnlik","is_probably unlikely",
                                    "is_probably unlikely"))))))))))))
# get rid of some columns
xdf.interpretation$comments<-xdf.interpretation$engagement<-xdf.interpretation$difficulty<-xdf.interpretation$language<-NULL
# add a value column, representing pairs <a,o> as a*10+0 (only exception is pair <10,10> represented as 1010)
xdf.interpretation$value=ifelse(xdf.interpretation$kind=="f",NA, #it makes sense only in guess-the-observation trials
                                ifelse(xdf.interpretation$access==10&xdf.interpretation$observation==10,1010,#takes care of exceptioncal case
                                       10*xdf.interpretation$access+xdf.interpretation$observation))
# save (to be used in analysis.r)
write.csv(xdf.interpretation,"clean_interpretation_data.csv", row.names = FALSE)

# visualize guess-the-state data (trials of kind "f"): discrete distributions over states
df.s=as.data.frame(xtabs(~expression+state,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="f"))))
# states as factor, it's really a categorical variable
df.s$state=as.factor(df.s$state)
colnames(df.s)=c("expression","state","count")
#plot
state.lines=ggplot(data=df.s)+geom_point(aes(x=state,y=count,group=expression),size=2.5)+geom_line(aes(x=state,y=count,group=expression))+facet_wrap(facets=~expression,ncol=3)+ scale_y_continuous(limits = c(0,70), breaks=seq(0,70,10))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
show(state.lines)

# visualize guess-the-observation 2dimensional data (trials of kind "c"): discrete joint distribution over pairs <a,o>
df.v=as.data.frame(xtabs(~value+expression,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))))
# fix type of value column
df.v$value=as.numeric(levels(df.v$value))[as.numeric(df.v$value)]
# add access and observation columns
df.v$access=as.factor(ifelse(df.v$value==1010,10,trunc(df.v$value/10)))
df.v$observation=as.factor(ifelse(df.v$value==1010,10,df.v$value-10*trunc(df.v$value/10)))
colnames(df.v)=c("value","expression","count","access","observation")
# colors
rf=colorRampPalette(rev(brewer.pal(11,'Spectral')))
r=rf(32)
# plots
value.heatmap = ggplot(df.v, aes(access,observation))+ geom_tile(aes(access,observation, fill = count))+ scale_fill_gradientn(colours=r)+facet_wrap(facets = ~ expression, ncol=3)+theme_bw()+theme(text=element_text( size=22), legend.position="right")
show(value.heatmap)

# visualize marginalized distributions over access and observation
# access
df.a=as.data.frame(xtabs(~expression+access,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))))
# accesss as factor, it's really a categorical variable
df.a$access=as.factor(df.a$access)
colnames(df.a)=c("expression","access","count")
# plot
access.lines=ggplot(data=df.a)+geom_point(aes(x=access,y=count,group=expression),size=2.5)+geom_line(aes(x=access,y=count,group=expression))+facet_wrap(facets=~expression,ncol=3)+ scale_y_continuous(limits = c(0,70), breaks=seq(0,70,10))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
# observation
df.o=as.data.frame(xtabs(~expression+observation,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))))
# observations as factor, it's really a categorical variable
df.o$observation=as.factor(df.o$observation)
colnames(df.o)=c("expression","observation","count")
# plot
observation.lines=ggplot(data=df.o)+geom_point(aes(x=observation,y=count,group=expression),size=2.5)+geom_line(aes(x=observation,y=count,group=expression))+facet_wrap(facets=~expression,ncol=3)+ scale_y_continuous(limits = c(0,70), breaks=seq(0,70,10))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
# show plots together
show(grid.arrange(access.lines,observation.lines, ncol = 2))


## compute bootstrapped 95% CIs for interpretation data
# state
states <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="f"))$state)) 
expressions <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="f"))$expression))  # experimental conditions
repetitions <- 1000 # how many samples?

# output is array expressions x state x repetitions, output1 is values x expressions x 3 (ci.low, mean, ci.high)
output <- array(dim=c(length(expressions),length(states),repetitions),dimnames = list(expressions,states,c(1:repetitions)))
output1 <- array(dim=c(length(expressions),length(states),3),dimnames = list(expressions,states,c("ci.low","mean","ci.high")))
for (e in 1:length(expressions)){ # for each condition...
  choices <- subset(subset(xdf.interpretation,xdf.interpretation$kind=="f"), subset(xdf.interpretation,xdf.interpretation$kind=="f")$expression==expressions[e])$state # get the vector with observed answers
  for (s in 1:length(states)){
    output[e,s,] <- bootstrap(x=choices,nboot=repetitions,theta = count,states[s])$thetastar
    output1[e,s,] <- c(ci.low(output[e,s,]),mean(output[e,s,]),ci.high(output[e,s,]))
  }
}
# as df
df.s.cis <- as.data.frame(cbind(melt(output1[,,1]),melt(output1[,,2])$value,melt(output1[,,3])$value))
colnames(df.s.cis) <- c("expression","state","ci.low","mean","ci.high")
# save
write.csv(df.s.cis,"complex_interpretation_state_cis.csv", row.names = FALSE)



# access
accesses <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))$access)) 
expressions <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))$expression))  # experimental conditions
repetitions <- 1000 # how many samples?

# output is array expressions x access x repetitions, output1 is values x expressions x 3 (ci.low, mean, ci.high)
output <- array(dim=c(length(expressions),length(accesses),repetitions),dimnames = list(expressions,accesses,c(1:repetitions)))
output1 <- array(dim=c(length(expressions),length(accesses),3),dimnames = list(expressions,accesses,c("ci.low","mean","ci.high")))
for (e in 1:length(expressions)){ # for each condition...
  choices <- subset(subset(xdf.interpretation,xdf.interpretation$kind=="c"), subset(xdf.interpretation,xdf.interpretation$kind=="c")$expression==expressions[e])$access # get the vector with observed answers
  for (a in 1:length(accesses)){
    output[e,a,] <- bootstrap(x=choices,nboot=repetitions,theta = count,accesses[a])$thetastar
    output1[e,a,] <- c(ci.low(output[e,a,]),mean(output[e,a,]),ci.high(output[e,a,]))
  }
}
# as df
df.a.cis <- as.data.frame(cbind(melt(output1[,,1]),melt(output1[,,2])$value,melt(output1[,,3])$value))
colnames(df.a.cis) <- c("expression","access","ci.low","mean","ci.high")
# save
write.csv(df.a.cis,"complex_interpretation_access_cis.csv", row.names = FALSE)


# observation
observs <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))$observation)) 
expressions <- levels(as.factor(droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))$expression))  # experimental conditions
repetitions <- 1000 # how many samples?

# output is array expressions x access x repetitions, output1 is values x expressions x 3 (ci.low, mean, ci.high)
output <- array(dim=c(length(expressions),length(observs),repetitions),dimnames = list(expressions,observs,c(1:repetitions)))
output1 <- array(dim=c(length(expressions),length(observs),3),dimnames = list(expressions,observs,c("ci.low","mean","ci.high")))
for (e in 1:length(expressions)){ # for each condition...
  choices <- subset(subset(xdf.interpretation,xdf.interpretation$kind=="c"), subset(xdf.interpretation,xdf.interpretation$kind=="c")$expression==expressions[e])$observation # get the vector with observed answers
  for (o in 1:length(observs)){
    output[e,o,] <- bootstrap(x=choices,nboot=repetitions,theta = count,observs[o])$thetastar
    output1[e,o,] <- c(ci.low(output[e,o,]),mean(output[e,o,]),ci.high(output[e,o,]))
  }
}
# as df
df.o.cis <- as.data.frame(cbind(melt(output1[,,1]),melt(output1[,,2])$value,melt(output1[,,3])$value))
colnames(df.o.cis) <- c("expression","observation","ci.low","mean","ci.high")
# save
write.csv(df.o.cis,"complex_interpretation_observation_cis.csv", row.names = FALSE)