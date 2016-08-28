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
comptemp<-read.table("tempcomp.txt",header=TRUE,sep='\t')
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


###############################################################################
#temporal competition
###############################################################################

#we start with binomial tests for each trial and we build a table of results
rez<-c()
for (i in c(1:12)) {
  temp<-binom.test(c(comptemp$Clone_1_BM[i],comptemp$Clone_2_BM[i]),p=0.5)
  rez<-rbind(rez,c(temp$parameter,temp$p.value))
  
}
rez


#we try to investigate if the effect of the competitive clone impact the 
#of the clone development
tempmod41<-glm(cbind(Clone_1_BM,Clone_2_BM)~Clone_2,family="binomial",
               data=comptemp[1:8,])
summary(tempmod41)
tempmod48<-glm(cbind(Clone_1_BM,Clone_2_BM)~Clone_2,family="binomial",
               data=comptemp[9:16,])
summary(tempmod48)
tempmod50<-glm(cbind(Clone_1_BM,Clone_2_BM)~Clone_2,family="binomial",
               data=comptemp[17:24,])
summary(tempmod50)


###############################################################################
#Metapopulation competition
###############################################################################

#very few things happened in the box, and we are mostly interested by the 
#dynamic within the pillbox. Therefore, we remove the "boite" lines
compspa<-compspa[compspa$Pilulier!="boite",]
compspa<-droplevels(compspa)

#let's check the evolution of the different type of aphid across time
op<-par(mfrow=c(3,1))
boxplot(compspa$L1L2~compspa$Day,main="Evolution temporelle Larve L1-L2",
        boxwex=0.5)
boxplot(compspa$L3L4~compspa$Day,main="Evolution temporelle Larve L3-L4",
        boxwex=0.5)
boxplot(compspa$Fem~compspa$Day,main="Evolution temporelle Femelle",
        boxwex=0.5)
par(op)

#we can also plot the evolution in for each experiment, the evoluation of each
#pillbox
op<-par(mfrow=c(3,3))

for (j in 1:9) {
  testid<-levels(compspa$Test)[j]
  limi<-max(compspa[compspa$Test==testid,c("L1L2","L3L4","Fem")])
  liad<-max(compspa[,("Fem")])
  plot(compspa$L1L2[compspa$Test==testid & 
                      compspa$Pilulier==levels(compspa$Pilulier)[1]]~
         compspa$Day[compspa$Test==testid & 
                       compspa$Pilulier==levels(compspa$Pilulier)[1]],type="l",
       col=1,ylim=c(0,limi),main=paste(testid,"- Evolution temporelle L1L2"),
       lwd=3,xlab="Nombre de jours",ylab="Effectif")
  for (i in 2:9) {
    lines(compspa$L1L2[compspa$Test==testid & 
                         compspa$Pilulier==levels(compspa$Pilulier)[i]]~
            compspa$Day[compspa$Test==testid & 
                          compspa$Pilulier==levels(compspa$Pilulier)[i]],type="l",
          col=i,ylim=c(0,limi),main=paste(testid,"- Evolution temporelle L1L2"),
          lwd=3,xlab="Nombre de jours",ylab="Effectif")
  }
  
  plot(compspa$L3L4[compspa$Test==testid & 
                      compspa$Pilulier==levels(compspa$Pilulier)[1]]~
         compspa$Day[compspa$Test==testid & 
                       compspa$Pilulier==levels(compspa$Pilulier)[1]],type="l",
       col=1,ylim=c(0,limi),main=paste(testid,"- Evolution temporelle L3L4"),
       lwd=3,xlab="Nombre de jours",ylab="Effectif")
  for (i in 2:9) {
    lines(compspa$L3L4[compspa$Test==testid & 
                         compspa$Pilulier==levels(compspa$Pilulier)[i]]~
            compspa$Day[compspa$Test==testid & 
                          compspa$Pilulier==levels(compspa$Pilulier)[i]],type="l",
          col=i,ylim=c(0,limi),main=paste(testid,"- Evolution temporelle L3L4"),
          lwd=3,xlab="Nombre de jours",ylab="Effectif")
  }
  
  plot(compspa$Fem[compspa$Test==testid & 
                     compspa$Pilulier==levels(compspa$Pilulier)[1]]~
         compspa$Day[compspa$Test==testid & 
                       compspa$Pilulier==levels(compspa$Pilulier)[1]],type="l",
       col=1,ylim=c(0,liad),main=paste(testid,"- Evolution temporelle Femelle"),
       lwd=3,xlab="Nombre de jours",ylab="Effectif")
  for (i in 2:9) {
    lines(compspa$Fem[compspa$Test==testid & 
                        compspa$Pilulier==levels(compspa$Pilulier)[i]]~
            compspa$Day[compspa$Test==testid & 
                          compspa$Pilulier==levels(compspa$Pilulier)[i]],type="l",
          col=i,ylim=c(0,liad),main=paste(testid,"- Evolution temporelle Femelle"),
          lwd=3,xlab="Nombre de jours",ylab="Effectif")
  }
}

par(op)


#let's check the distribution between pillbox at the end of the experiment
op<-par(mfrow=c(3,1))
boxplot(compspa$L1L2[compspa$Day==10]~compspa$Pilulier[compspa$Day==10],
        boxwex=0.5,main="Distribution par pilulier Larve L1-L2")
boxplot(compspa$L3L4[compspa$Day==10]~compspa$Pilulier[compspa$Day==10],
        boxwex=0.5,main="Distribution par pilulier Larve L3-L4")
boxplot(compspa$Fem[compspa$Day==10]~compspa$Pilulier[compspa$Day==10],
        boxwex=0.5,main="Distribution par pilulier Femelle")
par(op)


#plot of the final repartition of the clone in competition for
temp<-compspa[compspa$Day==10,]
levels(temp$Clone_1)<-c(1,2)
levels(temp$Clone_2)<-c(2,3)

op<-par(mfrow=c(3,3),mar=c(1,4.1,2,0.1))

for (i in 1:81) {
  if (!is.na(temp$BM_Clone_1)[i]) {
    pie(cbind(temp$BM_Clone_1,temp$BM.Clone_2)[i,],
        col=c(as.numeric(as.character(temp$Clone_1[i])),
              as.numeric(as.character(temp$Clone_2[i]))),
        main=temp$Unit[i],
        radius=sqrt(sum(temp[i,c("L1L2","L3L4","Fem")]))/20)
  } else {
    pie(5,col="transparent",main=temp$Unit[i],
        radius=sqrt(sum(temp[i,c("L1L2","L3L4","Fem")]))/20)
  }
  
}

par(op)


#model to explore if some factors affect the development of the clones
temp2<-cbind(temp,"totsum"=rowSums(temp[,c("L1L2","L3L4","Fem")]))
modefspa<-glm(totsum~Pilulier*Test,family="poisson",data=temp2)
summary(modefspa)
#no effect of the pillbox on the number of individual at the end of the 
#experiment (except for the pillbox 5 where the original individuals 
#were deposited). There is a lot of variability between Test

modefspa<-glm(totsum~Pilulier+Test+Clone_1+Clone_2,
              family="poisson",data=temp2)
summary(modefspa)

tapply(temp2$totsum,INDEX=temp2$Pilulier,FUN = mean)


