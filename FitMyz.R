###############################################################################
###############################################################################
#Script for the analysis of Fitness experiments
###############################################################################
###############################################################################

#first, setting the right working directory
setwd("~/Work/Rfichiers/Githuber/FitMyz_data")


###############################################################################
#Loading the dataset
###############################################################################

compspa<-read.table("compdata.txt",header=TRUE,sep='\t')
comptemp<-read.table("",header=TRUE,sep='\t')
single<-read.table("singledata.txt",header=TRUE,sep='\t')


###############################################################################
#Individual fitness traits
###############################################################################

#let see the structure of the data table
summary(single)
#and the distribution of the variable of interest
hist(single$age,by=single$Clone)

#we first build a complete model with all interaction included for the life 
#expectancy of the clone
#the explainatory variables are the Clone ID, the Posca treatment and whether 
#or not the larvae were removed along the experiment
modage<-glm(age~Clone*Posca*Larves_enlevees*date,family="poisson",data=single)
summary(modage)

#because the interactions do not show signigicant effect on the response 
#variable, we remove the interactions from the model
modage<-glm(age~Clone+Posca+Larves_enlevees+date,family="poisson",data=single)
summary(modage)
plot(modage)
modage<-glm(age~Clone+Posca+date,family="poisson",data=single)
summary(modage)
plot(modage)
modage<-glm(age~Clone+Posca,family="poisson",data=single)
summary(modage)
drop1(modage)
plot(modage)

library(lme4)
modrage<-glmer(age~Clone+Posca+Larves_enlevees+1|date,family="poisson",data=single)
modr<-lmer(age~Posca+Larves_enlevees+1|date,data=single)

single2<-single[single$date!="25-04-16",]
modage<-glm(age~Clone+Larves_enlevees,family="poisson",data=single2)
summary(modage)


modage<-glm(nb_larve~Clone+Posca+Larves_enlevees+date,family="gaussian",data=single2)
summary(modage)
drop1(modage,test="Chisq")
modage<-glm(nb_larve~Clone+Posca+Larves_enlevees,family="gaussian",data=single)
summary(modage)

modage<-lm(age~Clone,data=single)
summary(modage)
