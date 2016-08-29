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

compspatot<-read.table("compdata.txt",header=TRUE,sep='\t')
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
compspatot<-compspatot[compspatot$Pilulier!="boite",]
compspatot<-droplevels(compspatot)
#we remove the duplicate in the dataset
compspa<-compspatot[compspatot$dup=="y",]
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
    pie(cbind(temp$BM_Clone_1,temp$BM_Clone_2)[i,],
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


###############################################################################
#Metapopulation competition statistical analyses
###############################################################################

#with the number of individuals after 3 days
temp<-compspa[compspa$Day==3,]
#let sum up the total number of individuals, not taking into account 
#the stage
temp2<-cbind(temp,"totsum"=rowSums(temp[,c("L1L2","L3L4","Fem")]))

#let see if the pillbox number and the type of competition experiment 
#have an effect on the total number of aphids, the interaction are not 
#included because they are biologicaly irrelevant
modeftot<-glm(totsum~Pilulier+compID,family="poisson",data=temp2)
summary(modeftot)
modefL1L2<-glm(L1L2~Pilulier+compID,family="poisson",data=temp2)
summary(modefL1L2)
#After 3 days, we already detect difference between pillbox, but there is no
#clear pattern for which pillbox is more or less invaded. More importantly, 
#we already detect that competition with clone 50 show less individuals 
#(significant for one of the competition)

#with the number of individuals after 5 days
temp<-compspa[compspa$Day==5,]
#let sum up the total number of individuals, not taking into account 
#the stage
temp2<-cbind(temp,"totsum"=rowSums(temp[,c("L1L2","L3L4","Fem")]))

#let see if the pillbox number and the type of competition experiment 
#have an effect on the total number of aphids, the interaction are not 
#included because they are biologicaly irrelevant
modeftot<-glm(totsum~Pilulier+compID,family="poisson",data=temp2)
summary(modeftot)
modefL1L2<-glm(L1L2~Pilulier+compID,family="poisson",data=temp2)
summary(modefL1L2)
modefL3L4<-glm(L3L4~Pilulier+compID,family="poisson",data=temp2)
summary(modefL3L4)
#here we can see that both competition which include clone 50 are less 
#populated, even when you consider only the L3L4 larvae


#with the final number of individuals
temp<-compspa[compspa$Day==10,]
#let sum up the total number of individuals, not taking into account 
#the stage
temp2<-cbind(temp,"totsum"=rowSums(temp[,c("L1L2","L3L4","Fem")]))

#let see if the pillbox number and the type of competition experiment 
#have an effect on the total number of aphids, the interaction are not 
#included because they are biologicaly irrelevant
modefL1L2<-glm(L1L2~Pilulier+compID,family="poisson",data=temp2)
summary(modefL1L2)
modefL3L4<-glm(L3L4~Pilulier+compID,family="poisson",data=temp2)
summary(modefL3L4)
modefFem<-glm(Fem~Pilulier+compID,family="poisson",data=temp2)
summary(modefFem)
modeftot<-glm(totsum~Pilulier+compID,family="poisson",data=temp2)
summary(modeftot)
#at the end of the experiment, there are significantly less individuals 
#in the competition with clone 50 involved. It seems also that there are 
#a few more individuals in the pillbox located at the top of the box 
#(pillbox 1, 2, 3). Maybe related to some "environmental variation" in the 
#experimental setting. But there are a lot of variation among pillbox 
#between experiments

tapply(temp2$totsum,INDEX=temp2$Pilulier,FUN=mean)
tapply(temp2$totsum,INDEX=temp2$Pilulier,FUN=var)


#Let's check if the competition of the different clones affect the 
#final proportion of the clone. We include the pillbox in the model 
#because

compspatot10<-compspatot[compspatot$Day==10 & !is.na(compspatot$BM_Clone_1),]
compspatot10<-droplevels(compspatot10)
compspatot10<-cbind(compspatot10,
                    "totsum"=rowSums(compspatot10[,c("L1L2","L3L4","Fem")]))

spa41mod<-glmer(cbind(BM_Clone_1,BM_Clone_2)~Clone_2 + 1|Pilulier,
                family="binomial",
                data=compspatot10[compspatot10$Clone_1=='13 001 041',])
summary(spa41mod)

spa48mod<-glmer(cbind(BM_Clone_1,BM_Clone_2)~Clone_2 + 1|Pilulier,
                family="binomial",
                data=compspatot10[compspatot10$Clone_1=='13 001 048',])
summary(spa48mod)

spa50mod<-glmer(cbind(BM_Clone_1,BM_Clone_2)~Clone_2 + 1|Pilulier,
                family="binomial",
                data=compspatot10[compspatot10$Clone_1=='13 001 050',])
summary(spa50mod)

#it seems that there is no difference for the three clones when they are 
#competiting with the different clones




