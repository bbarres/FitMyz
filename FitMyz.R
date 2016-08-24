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
#we remove samples with no laid larvae or with clone dying during the first day
single<-single[single$age>1 & single$nb_larve>0,]

##life expectancy of the clone####
#and the distribution of the variable of interest
op<-par(mfrow=c(3,1))
plot(density(single$age[single$Clone=="13 001 041"]))
hist(single$age[single$Clone=="13 001 048"])
hist(single$age[single$Clone=="13 001 050"])
par(op)


#we first build a complete model with all interaction included for the life 
#expectancy of the clone
#the explainatory variables are the Clone ID, the Posca treatment and whether 
#or not the larvae were removed along the experiment
modage<-glm(age~Clone*Posca*Larves_enlevees,family="poisson",data=single)
summary(modage)

#we iteratively simplify the model by removing non-significant interactions
#and factors
modage<-glm(age~Clone+Posca+Larves_enlevees,family="poisson",data=single)
summary(modage)
#here is the final model
modagef<-glm(age~Clone+Posca,family="poisson",data=single)
summary(modagef)
plot(modagef)

##number of larvae produce
#we remove the repetition where the larvae weren't remove along the 
#experiment since they may have contributed to the total number of larvae
single2<-single[single$Larves_enlevees!="Non",]
#we first build a complete model with all interaction included for the 
#production of larvae
#the explainatory variables are the Clone ID, the Posca treatment and whether 
#or not the larvae were removed along the experiment
modnbl<-glm(nb_larve~Clone*Posca,family="poisson",data=single2)
summary(modnbl)

#we iteratively simplify the model by removing non-significant interactions
#and factors
modnbl<-glm(nb_larve~Clone+Posca+Larves_enlevees+Clone:Larves_enlevees+Clone:Posca,
            family="poisson",data=single)
summary(modnbl)
modnbl<-glm(nb_larve~Clone+Posca,family="poisson",data=single)
summary(modnbl)
modnbl<-glm(nb_larve~Clone+Posca+Larves_enlevees,family="poisson",data=single)
summary(modnbl)
#here is the final model
modnblf<-glm(nb_larve~Clone+Posca,family="poisson",data=single)
summary(modnblf)
plot(modnblf)



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
