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
library(bootstrap)
	

# source functions
source("useful_functions.r") # various useful definitions
source("prediction_functions.r") # the prediction functions, ie R implementation of the model 

# get data
df.production=read.csv("raw_production_data.csv",sep=";") # raw data from simple expression v2 production task
df.interpretation=read.csv("raw_interpretation_data.csv",sep=";") # raw data from simple expression interpretation task
df.observed.belief=read.csv("observed_belief.csv") # processed data from belief-about-the-urn task

## production task: clean and explore
xdf.production <- df.production

# how many subjects?
xdf.production$id %>% unique() %>% length()
# 89

# randomization of experimental conditions
xtabs(~value, xdf.production)
#value

# reported languages
xdf.production$language %>% levels()
#[1] "American English" "english"          "English"          "english "         "Russian"          "tamil"           
#[7] "Tamil"            "Turkish"          "Ukrainian"   

# drop non native english speaker
xdf.production=droplevels(subset(xdf.production,xdf.production$language=="American English" | 
                                   xdf.production$language=="English" |
                                   xdf.production$language=="english" |
                                   xdf.production$language=="english "))

# how many subjects left?
xdf.production$id %>% unique() %>% length()

# comments
xdf.production$comments %>% levels()
# interesting output:
# [7] "I think it has to be 100% to be certain "
# [8] "Interesting task! Love you guys. Mwah"
# [28] "To do certain things"

# reported difficulty
xtabs(~difficulty, xdf.production)/9 # 9 is the number of trials completed by each participant
# difficulty
# 0   1   2   3   4   5   6   7   8  # levels of difficulty  
# 29  18  14  9   1   8   3   1   1 

# produce and save a more compact, better looking data frame with less "useless" columns and a couple more useful columns
xdf.production$comments<-xdf.production$engagement<-xdf.production$difficulty<-xdf.production$language<-NULL
# add column with level of high-order uncertainty
xdf.production$uncertainty=ifelse(xdf.production$access==10,"none",ifelse(xdf.production$access>5,"low","high"))
# save (to be used in analysis.r)
write.csv(xdf.production,"clean_production_data.csv", row.names = FALSE)

# visualize expression choices in each condition with bar plots
# get counts of expression per condition
counts=t(xtabs(~answer+value,xdf.production))
# how many observation for each condition?
observations=as.vector(xtabs(~value,xdf.production))
# counts to percentage, as data frame
df.e=as.data.frame(melt(counts/observations*100))
colnames(df.e)=c("value","expression","percentage")
# add some columns to df for easier visualization
df.e$access=ifelse(df.e$value==1010,10,trunc(df.e$value/10))
df.e$observation=ifelse(df.e$value==1010,10,df.e$value-10*df.e$access)
# better looking labels for values
df.e$label=as.factor(paste0(df.e$observation," red balls out of ",df.e$access))
# better looking names for expressions
levels(df.e$expression)=c("certainly","certainly not","possibly","probably","probably not")
# change order of some factors to make plot more readable
df.e$label = factor(df.e$label,levels(df.e$label)[c(1,3,6,9,5,2,7,10,12,15,4,8,11,13,14)])
df.e$expression = factor(df.e$expression,levels(df.e$expression)[c(2,5,3,4,1)])
# plot
expression.bars=ggplot(data=df.e)+geom_bar(aes(x=expression,y=percentage,fill=expression),stat="identity")+
  scale_y_continuous(limits = c(0,100), breaks=seq(0,100,10))+coord_flip()+facet_wrap(facets = ~label, ncol=5)+
  theme_bw()+theme(strip.background = element_blank(),text=element_text(size=20), legend.position="none")
show(expression.bars)


## add bootstrapped 95% CIs
values <- levels(as.factor(xdf.production$value)) # experimental conditions
expressions <- levels(xdf.production$answer)
repetitions <- 1000 # how many samples?

# output is array values x expressions x repetitions, output1 is values x expressions x 3 (ci.low, mean, ci.high)
output <- array(dim=c(length(values),length(expressions),repetitions),dimnames = list(values,expressions,c(1:repetitions)))
output1 <- array(dim=c(length(values),length(expressions),3),dimnames = list(values,expressions,c("ci.low","mean","ci.high")))
for (v in 1:length(values)){ # for each condition...
    choices <- subset(xdf.production, xdf.production$value==values[v])$answer # get the vector with observed answers
    for (e in 1:length(expressions)){
      output[v,e,] <- bootstrap(x=choices,nboot=repetitions,theta = percentage,expressions[e])$thetastar
      output1[v,e,] <- c(ci.low(output[v,e,]),mean(output[v,e,]),ci.high(output[v,e,]))
  }
}
# as df
df.e.cis <- as.data.frame(cbind(melt(output1[,,1]),melt(output1[,,2])$value,melt(output1[,,3])$value))
colnames(df.e.cis) <- c("value","message","ci.low","mean","ci.high")
# save
write.csv(df.e.cis,"simple_production_cis.csv", row.names = FALSE)


## Classic statistics
# multinomial logistic regression on answer choices
# mlogit needs data in long form
tempdata=xdf.production
# rearrange, so that possibly is reference (mlogit.data seems to force alphabetic order)
levels(tempdata$answer)=c("E_certainly","B_certainlyNot","A_possibly","D_probably","C_probablyNot")
long.tempdata0 = mlogit.data(tempdata,shape="wide",choice="answer")
m0 = mlogit(answer ~ 0|observation+access, data=long.tempdata0)
summary(m0)
# wrt possibly, access and observation have significant effects on every message

# does including interaction between access and observation result in a better model?
m0i = mlogit(answer ~ 0|observation*access, data=long.tempdata0)
# model comparison
AIC(m0,m0i)
# df      AIC
# m0  12 1591.128
# m0i 16 1593.425
# nope, barely distinguishable

# similarly to what we did in cogsci2016 paper
# we can AIC-compare this model to a model with only observed proportion as predictor...
tempdata1=tempdata
tempdata1$proportion=tempdata1$observation/tempdata1$access
long.tempdata1 = mlogit.data(tempdata1,shape="wide",choice="answer")
m1 = mlogit(answer ~ 0|proportion, data=long.tempdata1)
summary(m1)
# wrt possibly, proportion has significant effect on every message
# model comparison
AIC(m0,m1)
# df      AIC
# m0 12 1591.128
# m1  8 1729.504
# --> despite added complexity, m0 is AIC-better than m1

# ...or to a model with only point estimates (eg mode or EV) of objective chance obtained by measured beliefs
# first, take subset of observed beliefs df picking only the values <a,o> that we care about
xdf.observed.belief=droplevels(subset(df.observed.belief,df.observed.belief$value %in% xdf.production$value))
# then add the point estimates to each row of the full df, as function of corresponding value
for(v in 1:length(xdf.production[,1])){
  xdf.production$obs_bel_mode[v]=c(0:10)[which.max(subset(xdf.observed.belief,xdf.observed.belief$value==xdf.production$value[v])$mean)]
  xdf.production$obs_bel_ev[v]=sum(c(0:10)*subset(xdf.observed.belief,xdf.observed.belief$value==xdf.production$value[v])$mean)
}
# compute the two models: explain answer choice given mode of belief distribution or expected value
tempdata2=xdf.production
long.tempdata2 = mlogit.data(tempdata2,shape="wide",choice="answer")
m2 = mlogit(answer ~ 0|obs_bel_mode, data=long.tempdata2)
summary(m2)
m3 = mlogit(answer ~ 0|obs_bel_ev, data=long.tempdata2)
summary(m3)
# model comparison
AIC(m0,m2,m3)
# df      AIC
# m0 12 1591.128
# m2  8 1689.593
# m3  8 1620.963
# --> m0 always the AIC-best


## MLE to compute AIC given data
## basic set up of the model
n=10 # balls in the urn
states=c(0:n) # state space, ie the possible numbers of red balls in the urn
# full set of possible values, ie natural numers coding <a,o> pairs, eg <8,3> coded as 83
# generated calling function from useful_function.r
V=valuesF()$pairValues
A=valuesF()$accessValues
O=valuesF()$observationValues
access=unique(A)
observation=unique(O)
# values used in the production exp, ie experimental conditions
V.prod=df.production$value %>% unique() %>% sort() # pairs <a,o>
A.prod=ifelse(V.prod==1010,10,trunc(V.prod/10)) # access values
O.prod=ifelse(V.prod==1010,10,V.prod-(A.prod*10)) # observation values
if (V.prod[length(V.prod)]==110){V.prod[length(V.prod)]=1010}
access.prod = unique(A.prod)
observation.prod = unique(O.prod)
# messages
messages=levels(df.production$answer)

# massage data into right form
d = melt(with(xdf.production, table(answer, access, observation)))
# generate predictions with dummy values
predDummy = predictions.sp(1,1,1,0.5,0.5,0.5)$speaker
predDummy$access=ifelse(predDummy$value==1010,10,trunc(predDummy$value/10))
predDummy$observation=ifelse(predDummy$value==1010,10,(predDummy$value-(predDummy$access*10)))
# extract relevant factors from predictions
access = levels(factor(predDummy$access))
observation = levels(factor(predDummy$observation))
answer = levels(factor(predDummy$answer))
obs = 0
total = 0

for (i in 1:nrow(predDummy)) {
  dSub = filter(d, answer == predDummy$answer[i], 
                access == predDummy$access[i],
                observation == predDummy$observation[i])
  obs[i] = dSub$value
  dSub = filter(d, observation == predDummy$observation[i], 
                access == predDummy$access[i])
  total[i] = sum(dSub$value)
}

nLLOptim = function(par) {
  alpha = par[1]
  beta = par[2]
  lambda = par[3]
  theta.possibly = par[4]
  theta.probably = par[5]
  theta.certainly= par[6]
  if (alpha < 1 | beta < 1 | theta.probably <= 0.5 | theta.possibly <0 | theta.possibly >= 1 | theta.probably >= 1 | theta.certainly >= 1){
    NA
  } else {
    pred.sp = predictions.sp(alpha,beta,lambda,theta.possibly,theta.probably,theta.certainly)$speaker
    pred=pred.sp$probability
    - sum(dbinom(obs, total, pred, log = TRUE)) 
  }
}
fitOpt = optim(par = c(1,1,5,0.55,0.55,0.55), nLLOptim)

show(paste("Best-fit parameters: ", 
           paste(c("alpha =","beta =","lambda =","theta.possibly =","theta.probably =","theta.certainly ="), round(fitOpt$par, 3), sep = " ", collapse = " ")))
show(paste("AIC of best fit: ", round( fitOpt$value*2 + 8, 3)))
# [1] "Best-fit parameters:  alpha = 1.49 beta = 1.994 lambda = 5.508 theta.possibly = 0.201 theta.probably = 0.601 theta.certainly = 0.933"
# [1] "AIC of best fit:  649.931"



############################################################################


## interpretation task: clean and explore
xdf.interpretation <- df.interpretation

# how many subjects?
xdf.interpretation$id %>% unique() %>% length()
# 147

# reported languages
xdf.interpretation$language %>% levels()
# [1] "A1F9KLZGHE9DTA" "eng"            "Engish"         "Englisg"        "english"        "English"        "ENGLISH"       
# [8] "English "       "englisj"        "Turkish" 


# drop non native english speakers
xdf.interpretation=droplevels(subset(xdf.interpretation,xdf.interpretation$language!="A1F9KLZGHE9DTA" & xdf.interpretation$language!="Turkish"))

# how many subjects left?
xdf.interpretation$id %>% unique() %>% length()
# 145

# comments
xdf.interpretation$comments %>% levels()
# interesting output:
# [6] "I had to stop and think about all the possibilities. I have most likely overestimated in some instances, but my responses should still fit the bill. Thanks."  
# [8] "In the future, please put the overall scale back on the remaining screens after the warm-up. I'm fairly certain that the 'probable' as the center and 'possible' was just below it, but I'm not COMPLETELY certain. This is important because people will estimate in ranges of 20% (there were 5 choices) if they are rational. Or even if they are not, numbers and probabilities still dictate our intuitive responses, contra aesthetic theorists. Good luck with the results!"
# [10] "Interesting but completely random task unless more than 5 balls are taken all with the same color " 
# [15] "It was hard guessing like that. It really made me think. Interesting task, though."  
# [41] "That was interesting because there were so many ways to look at it." 
# [43] "The task was pretty confusing." 
# [46] "This was interesting, but confusing for me, too."  

# reported difficulty
xtabs(~difficulty, xdf.interpretation)/10 # 10 is the number of trials completed by each participant
# difficulty
#0  1  2  3  4  5  6  7  8  9 10 
#28  7 20 26  7 27  6  9  8  3  3 
#ok

# has anybody said acc=0 and obs=0 for any expression?
xdf.interpretation %>% filter(access==0 & observation==0) %>% select(id) %>% n_distinct()
# 1
# drop it
pessimists=levels(droplevels(subset(xdf.interpretation, xdf.interpretation$access==0 & xdf.interpretation$observation==0))$id)
xdf.interpretation=droplevels(subset(xdf.interpretation,!xdf.interpretation$id %in% pessimists))
# how many subjects left?
xdf.interpretation$id %>% unique() %>% length()
# 144

# produce and save a more compact, better looking data frame
# rename expression levels, to make them more readable and match names in production experiment and model
xdf.interpretation$expression=ifelse(xdf.interpretation$expression=="cer","certainly",
                       ifelse(xdf.interpretation$expression=="cerNot","certainlyNot",
                              ifelse(xdf.interpretation$expression=="poss","possibly",
                                     ifelse(xdf.interpretation$expression=="prob","probably","probablyNot"))))

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
max=max(df.s$count)
state.lines=ggplot(data=df.s)+geom_point(aes(x=state,y=count,group=expression),size=2.5)+geom_line(aes(x=state,y=count,group=expression))+facet_wrap(facets=~expression,ncol=5)+ scale_y_continuous(limits = c(0,max), breaks=seq(0,max,15))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
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
value.heatmap = ggplot(df.v, aes(access,observation))+ geom_tile(aes(access,observation, fill = count))+ scale_fill_gradientn(colours=r)+facet_wrap(facets = ~ expression, ncol=5)+theme_bw()+theme(text=element_text( size=22), legend.position="right")
show(value.heatmap)

# visualize marginalized distributions over access and observation
# access
df.a=as.data.frame(xtabs(~expression+access,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))))
# accesss as factor, it's really a categorical variable
df.a$access=as.factor(df.a$access)
colnames(df.a)=c("expression","access","count")
# plot
access.lines=ggplot(data=df.a)+geom_point(aes(x=access,y=count,group=expression),size=2.5)+geom_line(aes(x=access,y=count,group=expression))+facet_wrap(facets=~expression,ncol=5)+ scale_y_continuous(limits = c(0,70), breaks=seq(0,70,10))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
# observation
df.o=as.data.frame(xtabs(~expression+observation,data=droplevels(subset(xdf.interpretation,xdf.interpretation$kind=="c"))))
# observations as factor, it's really a categorical variable
df.o$observation=as.factor(df.o$observation)
colnames(df.o)=c("expression","observation","count")
# plot
observation.lines=ggplot(data=df.o)+geom_point(aes(x=observation,y=count,group=expression),size=2.5)+geom_line(aes(x=observation,y=count,group=expression))+facet_wrap(facets=~expression,ncol=5)+ scale_y_continuous(limits = c(0,70), breaks=seq(0,70,10))+theme_bw() +theme(strip.background = element_blank(),text=element_text( size=24), legend.position="bottom")
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
write.csv(df.s.cis,"simple_interpretation_state_cis.csv", row.names = FALSE)



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
write.csv(df.a.cis,"simple_interpretation_access_cis.csv", row.names = FALSE)


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
write.csv(df.o.cis,"simple_interpretation_observation_cis.csv", row.names = FALSE)

