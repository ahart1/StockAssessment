#Notes
#standardization of SSB and Y (/1000) was done to convert data to same units, may not be necessary always

#This is the Function that calculates Abundance and Catch
AgeBasedMSY <- function(data){

  #List of Variables at age and starting parameters
  Age<-data$Age            #Age 1-20
  F2007<-data$F.2007       #Fishing mortality 2007 bluefin tuna data
  Ma<-data$Maturity.ma     #Mortality-at-age
  Pa<-data$Pa              #Partial recruitment-at-age
  La<-data$Length.La       #Length-at-age
  Wa<-data$Weight.wa       #Weight-at-age
  k<-data$k                #growth parameter
  Linf<-data$Linf          #maximum length growth curve parameter
  t0<-data$t0              #initial size growth curve parameter
  Ffull<-data$Ffull        #Full fishing mortality 
  alpha<-data$alpha        
  beta<-data$beta
  N1<-data$N1              #Starting abundance of age 1 fish in year 1 of model run
  F<-data$F                #Fishing mortality (may vary for a range of results)
  b1<-data$b1              #Parameter for recruitment curve
  b2<-data$b2              #Parameter for recruitment curve
  M<-data$M                #Natural Mortality
  
  
  ########Abundance Calculations Begin Here##############
  
  #Abundance for Year 1
  #Set up matrix 50 row(year) by 20 column(age)
  N <- matrix(rep(NA,50*20),50,20)
  #Compute abundance year 1
  N[1,1]<-N1[1]
  for(a in 1:19){
    N[1,a+1]<-N[1,a]*exp(-(Pa[a]*F[1]+M[1]))
  }
  
  #Create list for SSB and R to fill in with calculated values below
  SSB<-rep(0,50)
  R<-rep(0,49)
  
  
  #Abundance for year 2-50
  #This is code for year 2-50, it begins by calculating SSB for year 1 and recruitment which is used as N1 for year 2
  for(t in 1:50){
    #Define SSB for each year(t) and age(a), standardize /1000
    total<-0
    for(a in 1:20){
      total <- total+ N[t,a]*Ma[a]*Wa[a]
    }
    SSB[t]<-total/1000
    
    #Stop loop for recruitment calculations in last year(year 50) so no abundance for age 1 in year 51 calculated
    if(t==50) {break}
    
    #Recruitment(R)
    R[t]<-b1[1]*SSB[t]/(SSB[t]+b2[1])
    
    #Set age 1 abundance equal to recruits in all but year 1
    N[t+1,1]<-R[t]
    
    #Calculate Abundance-at-age for year 2-50
    for(a in 1:19){
      N[t+1,a+1]<-N[t,a]*exp(-(Pa[a]*F[1]+M[1]))
    }
  }
  
  print(N)
  print(R)
  print(SSB)
  
  ############Catch at Age Calculations Begin Here###########
  
  #Set up matrix 50 row(year) by 20 column(age)
  C <- matrix(rep(NA, 50*20), 50, 20)
  
  #Set up list for Yield values
  Y <- rep(0,50)
  
  
  
  for(t in 1:50){
    #Make temporary variable (total) for use in yield calculations
    total <- 0
    
    #Calculate catch-at-age
    for(a in 1:20){
      C[t,a] <- N[t,a]*Pa[a]*F[1]*((1-exp(-(M[1]+Pa[a]*F[1])))/(M[1]+Pa[a]*F[1]))
      
      #Yield Calculations (sum of products from age 1 to 20, standardize by dividing by 1000)
      total <- total+C[t,a]*Wa[a]
      
    }
    Y[t] <- total/1000
  }
  
  print(C)
  print(Y)

  return(list(N=N,R=R,SSB=SSB,C=C,Y=Y))
}






##############Use Function to Obtain Abundance, Recruitment, Catch, Yield##########################
setwd("/Users/arhart/Documents")

#Read in Data File
InitialData<-read.csv("ModelPracticeData.csv", header=TRUE)

#Read return info into Output
Output <- AgeBasedMSY(InitialData)

print(Output$N)
print(Output$R)
print(Output$SSB)
print(Output$C)
print(Output$Y)







############Calculate Equilibrium SSB, Yield, Recruitment (Avg over yr 26-50)##################

#Equilibrium Recruitment
#R was set equal to N[t+1,1] (abundance of age 1 in the following year from which SSB came) there are values of N[t+1,1] for year 26 to 50 so the average of these values is equivalent to equilibrium recruitment
EquilibriumR <- mean(Output$N[26:50,1])
print(EquilibriumR)

#Equilibrium Yield
EquilibriumY <- mean(Output$Y[26:50])
print(EquilibriumY)

#Equilibrium Spawning Stock Biomass
EquilibriumSSB <- mean(Output$SSB[26:50])
print(EquilibriumSSB)

