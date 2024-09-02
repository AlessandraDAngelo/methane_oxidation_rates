#R script for calculating the microbial methane oxidation rates of in-vitro samples, following: 1st orded kinetics, Mahieu et al. (2006), Uhlig and Loose (2017). 
#D'Angelo A. and Loose B. (2020) - University of Rhode Island, Graduate School of Oceanography

#load dataframe "data"

library(dplyr) 
library(ggplot2)

#subset, based on samples
DF = data.frame()

Ubags = unique(data$bag.Num)

for (B in Ubags){
  
  BagID = data$bag.Num == B
  
  BagSub = data[BagID,]
  
  #calculate the methane concentration in the water
  a=1027 #g/L
  Vw=(BagSub$bag.M/a) #L
  BagSub$Conc_w<- BagSub$tot_mass/Vw #mol/L
  
  # Calculate the methane oxidation rate constant for concentration
  BagSub$ln.cCH4.water<- log(BagSub$Conc_w)
  BagSub$ln_tot.mass<- log(BagSub$tot_mass)
  # Calculate the linear regression for k from methane concentration
  BagSub$k_ox=lm(BagSub$ln_tot.mass ~ BagSub$Date)$coef[2]
  BagSub$k_ox= BagSub$k_ox*-1 #time^-1

  # Calculate the oxidation rates from isotopic fractionation
  #give an index for the run_num again
  Run<- unique(BagSub$run_num)
  #set Run as a list
  r<- list(Run)
 
  # From Mahieu et al. (2006) and the Rayleigh model: ln(CH4/CH4.t0)= a/1-a ln(1000+d13C/1000+d13C.t0), where a can be obtained by linear regression of experimental data.
  
  for (j in r) {
    ln.B= log((1000+BagSub$iso.data.calibrated)/(1000+BagSub$iso.data.calibrated[1]))
    
    # Fractionation factor alpha=k12/k13
    
    # From Uhlig and Loose (2017)
    alpha.upper.bound=1.025
    # From our calculation below
    alpha.lower.bound=1.008 
    
    reg.alpha <- lm(ln.B ~BagSub$ln.cCH4.water)
    lm.alpha <- summary(reg.alpha) # details on linear regression
    m.alpha <- lm.alpha$coefficients[2,1]
    term.alpha.up=  alpha.up/(1-alpha.up) 
    term.alpha.low= alpha.up/(1-alpha.low) 
    
    BagSub$term.left.up=ln.B*term.alpha.up
    BagSub$term.left.low=ln.B*term.alpha.low
    
    #lower bound
    BagSub$kox.delta.low= lm(BagSub$term.left.low ~ BagSub$Date)$coef[2] #time^-1
    #negative slope
    BagSub$kox.delta.low<-  BagSub$kox.delta.low*-1
    
    #upper bound
    BagSub$kox.delta.up= lm(BagSub$term.left.up ~ BagSub$Date)$coef[2] #time^-1
    #negative slope
    BagSub$kox.delta.up<-  BagSub$kox.delta.up*-1
  }
  DF <- rbind(DF,BagSub)
}

#create a quality flag for the kox statistically significant
data$QC_kox.a<- ifelse(data$k_ox>0 & data$kox.delta.up<0,0, 1)
data$QC_kox.b<-  ifelse(data$k_ox<0 & data$kox.delta.low>0,0, 1)

data.QC<- subset(data, QC_kox.a == 1 & QC_kox.b == 1)

# Caculate the methane microbial rate of oxidation by adding the column with CH4 concentrations from in-situ samples. 
data$r_ox = data$k_ox * data$CH4.conc #mol/L d 

# Plot kox vs kox_delta 
ggplot(NULL)+ 
  geom_point(data = MOSAiC.methox.results.QC, aes(x= k_ox, y =kox.delta.up), pch = 21, fill = "blue", color="dark blue", size = 5)+
  geom_smooth(data = MOSAiC.methox.results.QC, aes(x= k_ox, y =kox.delta.up), method = "lm", color = "red", se = FALSE, size = 1)+
  geom_abline(color="red",linetype="dashed", size=1)+
  theme_bw()+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(hjust = .5, size = 15), axis.text.y = element_text(size = 15)) +
  labs (title = "Methane oxidation rate constants")+
  xlab (expression(paste("k" ["ox.mass.balance"], (d^"-1"))))+
  ylab (expression(paste("k" ["ox.isotope.ratio"],(d^"-1"))))+
  theme(legend.position =  "right") 
library(dplyr)

#References:
# Mahieu, K., A. D. Visscher, P. A. Vanrolleghem, and O. V. Cleemput. 2006. Carbon and hydrogen isotope fractionation by microbial methane oxidation: Improved determination. Waste Manag. 26: 389–398. doi:10.1016/j.wasman.2005.11.006.
# Uhlig, C., Loose, B. 2017. Using stable isotopes and gas concentrations for independent constraints on microbial methane oxidation at Arctic Ocean temperatures. Limnol. Oceanogr. Methods, 15, 737-751. https://doi.org/10.1002/lom3.10199.
