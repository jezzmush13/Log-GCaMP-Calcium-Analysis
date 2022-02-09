# Log-GCaMP-Calcium-Analysis

#GCaMP Calcium Analysis log normal

#GCaMP Calcium Analysis JMR

#load the library 
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)

#load the data 
filenames <- list.files(pattern="*.csv") #read all the .csv files 


#create a dataframe; Remove missing values
df<-purrr::map_df(filenames, read_csv, .id = 'filename') 
df<- na.omit(df) #remove missing values 

#change the names to lowercase
df$roiName<- tolower(df$roiName)
df$Spot<- tolower(df$Spot)
df$Condition<- tolower(df$Condition)
df$drug<- tolower(df$drug)

#Change nominal and date format
df$filename = factor(df$filename)
df$roiName = factor(df$roiName)
df$Animal = factor(df$Animal)
df$drug = factor(df$drug)
df$Spot = factor(df$Spot)
df$filename = factor(df$filename)
df$date = factor(df$date)
view(df)

#Modifying variables and cleaning data frame
#Short df and create a unique ROIname
df_short <- df[ ,c('date', 'drug', 'Condition', 'Spot', 'Animal', 'roiName', 'amplitude', 'halfWidth','peakAUC')]
#Replace after drugs for baseline
df_short$drug <- str_replace_all(df_short$drug, 'after drugs', 'baseline')
#adding unique Spot
df_short$Unique_Spot <-paste(df_short$Animal, df_short$Spot, sep= "_")
#adding unique ROI name
df_short$Unique_ROIname <-paste(df_short$Animal, df_short$Spot, df_short$roiName, sep= "_")
#adding Duration of event
df_short$Duration <-df_short$halfWidth*2
#Replacing the labiling of "link" to process in pericyte structure
df_short$roiName <- str_replace_all(df_short$roiName, 'l1', 'plink')
#create a new variable for process and soma
df_short$ROI<- ifelse(grepl(pattern = "p", df_short$roiName, ignore.case = T),"Process","Soma")
#verifying categorical variables of the data frame
unique(df_short$drug)
unique(df_short$ROI)
unique(df_short$Animal)
df_short <- df_short[ ,c('date', 'Animal', 'Unique_Spot', 'Unique_ROIname', 'ROI', 'Condition', 'drug', 'amplitude', 'Duration','peakAUC')]
view(df_short)

#Remove negative amplitudes before calculating the mean
df_short<- df_short%>%
  filter(amplitude>0)


#Mean of each trial by ROI 

# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(df_short, c("date", 'Animal', "Unique_Spot", "Unique_ROIname", "ROI", "Condition", "drug"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(Duration,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude))


view(ROI.means)

#Converting to nominal factor the modified data frame
ROI.means$date = factor(ROI.means$date)
ROI.means$Animal = factor(ROI.means$Animal)
ROI.means$Unique_Spot= factor(ROI.means$Unique_Spot)
ROI.means$Unique_ROIname = factor(ROI.means$Unique_ROIname)
ROI.means$ROI = factor(ROI.means$ROI)
ROI.means$Condition = factor(ROI.means$Condition)
ROI.means$drug = factor(ROI.means$drug)



########################################################################################################
#STATISTICS

#STATISTICS

summary(ROI.means) 

#Amplitude

# explore the amplitude
library(plyr)
ddply(ROI.means, ~ drug * Condition, function(data) summary(data$amp_mean))
ddply(ROI.means, ~ drug * Condition, summarise, amp.mean=mean(amp_mean), amp.sd=sd(amp_mean))

#MODIFICATION OF THE DATA WITH LOGNORMAL DISTRIBUTION
#creating log of amp_mean
ROI.means$logamp_mean = log(ROI.means$amp_mean) # log transform
View(ROI.means) # verify

hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "stim",]$logamp_mean)

boxplot(logamp_mean ~ drug * Condition, data=ROI.means, xlab="Drug.Condition", ylab="amplitude") # boxplots
with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction

#Testing for normality in the residuals
library(plyr)
m = aov(logamp_mean ~ drug * Condition, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(logamp_mean ~ drug * Condition, data=ROI.means, center=mean) # Levene's test
leveneTest(logamp_mean ~ drug * Condition, data=ROI.means, center=median) # Brown-Forsythe test


# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"


#Anova of the different LMM on amplitude One categorical factor, with unique ROI name as random effect.
m1 <- lmer(logamp_mean ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
m2 <- lmer(logamp_mean ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
m3 <- lmer (logamp_mean ~ (ROI )  + (1|Unique_ROIname), data=ROI.means)
m4 <- lmer(logamp_mean ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
m5 <- lmer(logamp_mean ~ (drug * ROI)  + (1|Unique_ROIname), data=ROI.means)
m6 <- lmer(logamp_mean ~ (Condition * ROI)  + (1|Unique_ROIname), data=ROI.means)
amp_anova <- anova(m1, m2, m3, m4, m5, m6, refit = FALSE)
print(amp_anova)


#LMM with Tukey's correction. 
m4 = lmer(logamp_mean ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
emmeans(m, pairwise ~ drug*Condition)
#other option
m7 <- lsmeans(m, pairwise ~ drug*Condition, adjust="tukey")
summary(m7)


#LMM with holm-sequential Bonferroni 
# LMM with Spot as random effect for DRUG factor
m = lmer(logamp_mean ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for Condition factor
m = lmer(logamp_mean ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ Condition)), test=adjusted(type="holm"))


# LMM with Spot as random effect drug*Condition
m = lmer(logamp_mean ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?





