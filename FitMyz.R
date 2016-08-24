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
#distribution of the variable of interest
hist(single$age)
op<-par(mfrow=c(3,1))
hist(single$age[single$Clone=="13 001 041"])
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

barplot(c(mean(single$age[single$Clone=="13 001 041"]),
          mean(single$age[single$Clone=="13 001 048"]),
          mean(single$age[single$Clone=="13 001 050"])))


##number of larvae produce####
#we remove the repetition where the larvae weren't remove along the 
#experiment since they may have contributed to the total number of larvae
single2<-single[single$Larves_enlevees!="Non",]

#distribution of the studied variable
hist(single2$nb_larve)
op<-par(mfrow=c(3,1))
hist(single2$age[single2$Clone=="13 001 041"])
hist(single2$age[single2$Clone=="13 001 048"])
hist(single2$age[single2$Clone=="13 001 050"])
par(op)

#we first build a complete model with interaction included for the 
#production of larvae
#the explainatory variables are the Clone ID and the Posca treatment
modnbl<-glm(nb_larve~Clone*Posca,family="poisson",data=single2)
summary(modnbl)
plot(modnbl)

#modnbl<-glm(nb_larve~Clone+Posca,family="poisson",data=single2)

barplot(c(mean(single2$nb_larve[single2$Clone=="13 001 041"]),
        mean(single2$nb_larve[single2$Clone=="13 001 048"]),
        mean(single2$nb_larve[single2$Clone=="13 001 050"])))
