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

#we first build a complete model with all interaction included for the life 
#expectancy of the clone
#the explainatory variables are the Clone ID, the Posca treatment and wether 
#or not the larvae were removed along the experiment
modage<-glm(nb_larve~Clone*Posca*Larves_enlevees*date,family="gaussian",
            data=single)
summary(modage)

#because the interactions do not show signigicant effect on the response 
#variable, we remove the interactions from the model
modage<-glm(nb_larve~Clone+Posca+Larves_enlevees+date,family="gaussian",data=single)
summary(modage)


