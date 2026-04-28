#iSoM Analysis Script 

library(plyr)
library(Rmisc)
library(ggplot2)
library(dplyr)
library(mediation)
library(JSmediation)
library(ggsignif)
library(tidyverse)
library(stringr)
library(wrapr)
library(reshape2)
library(purrr)
library(Hmisc)
library(corrplot)
library(ggpubr)
library(interactions)

#d' function
d.prime <- function(x) {
  x$hit.N <- x$hit + x$miss
  x$FA.N <- x$FA + x$CR
  x$hit.rate <- x$hit/x$hit.N
  x$FA.rate <- x$FA/x$FA.N
  x$hit.rate <- ifelse(x$hit.rate == 0, (1/(2 * (x$FA.N + x$hit.N))), x$hit.rate)
  x$hit.rate <- ifelse(x$hit.rate == 1, 1 - (1/(2 * (x$FA.N + x$hit.N))), x$hit.rate)
  x$FA.rate <- ifelse(x$FA.rate == 0, (1/(2 * (x$FA.N + x$hit.N))), x$FA.rate)
  x$FA.rate <- ifelse(x$FA.rate == 1, 1 - (1/(2 * (x$FA.N + x$hit.N))), x$FA.rate)
  x$zhit.rate <- qnorm(x$hit.rate)
  x$zFA.rate <- qnorm(x$FA.rate)
  x$d.prime <- x$zhit.rate -  x$zFA.rate 
  x$bias <- (x$zhit.rate + x$zFA.rate)/2
  x$A.prime <- (0.5 + ((x$hit.rate-x$FA.rate)*(1 + x$hit.rate - x$FA.rate))/(4*x$hit.rate*(1-x$FA.rate)))
  x$A.normalized <- 2*asin(sqrt(x$A.prime))
  return(x)
}

###HBD Task Scoring###

#Read in HBD Task Data
#Generate List of Files to Analyze
hbd_files <- list.files(path = "/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Behavioral",
                        pattern = ".csv",
                        full.names = TRUE)

#Read in all the the files and add column that lists file name 
hbd.raw <- hbd_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = TRUE), .id = "participant_id")

#Clean up participant names
hbd.raw$participant_id <- str_remove(hbd.raw$participant_id, "/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Behavioral/")
hbd.raw$participant_id <- str_remove(hbd.raw$participant_id, "_HBD_Output.csv")

#Label Response Types (i.e., Hit, Miss, FA, CR)
hbd.raw <- within(hbd.raw,{
  response.type <- ifelse(Delay=="200" & Response == "1", "hit",
                          ifelse(Delay=="200" & Response =="0", "miss",
                                 ifelse(Delay=="500" & Response=="1","FA",
                                        ifelse(Delay=="500" & Response=="0", "CR", "NA"))))
})

#Tally Responses per Participant
hbd.tally <- dcast(data=hbd.raw, participant_id ~ response.type, fun.aggregate = length, value.var = "response.type", fill = 0)

#Calculate d prime and bias (c)
hbd.dprime1 <- d.prime(hbd.tally)

#Add Mean Confidence
hbd.conf <- summarySE(data = hbd.raw, measurevar = "Confidence", groupvars = "participant_id")
hbd.conf <- hbd.conf[c("participant_id", "Confidence")]

hbd.dprime <- merge(hbd.dprime1, hbd.conf, by = "participant_id")

#Fix isom_p05 id
hbd.dprime$participant_id <- ifelse(hbd.dprime$participant_id == "isom_p05b", "isom_p05", hbd.dprime$participant_id)
hbd.dprime$participant_id <- ifelse(hbd.dprime$participant_id == "isom_p44b", "isom_p44", hbd.dprime$participant_id)
hbd.dprime$participant_id <- ifelse(hbd.dprime$participant_id == "isom_p32b", "isom_p32", hbd.dprime$participant_id)

#Calculate good detectors
# --- Extended Binomial function  --------------------
hbd.dprime$correct <- hbd.dprime$hit + hbd.dprime$CR
hbd.dprime$accuracy <- hbd.dprime$correct/40

#hbd.dprime$p_value <- binom.test(hbd.dprime$correct, 40, p = 0.5, alternative = "greater")$p.value
hbd.dprime$p_value <- mapply(function(x) {
  binom.test(x, 40, p = 0.5, alternative = "greater")$p.value
}, hbd.dprime$correct)

hbd.dprime$good_detector <- hbd.dprime$accuracy >= 0.6 & hbd.dprime$p_value < 0.05



###Physio Data###
isom_hr <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/isom_mean_hr.csv")
isom_hr$participant_id <- str_remove(isom_hr$filename, "_physio.acq")
isom_hr$meanIBI <- 60000/isom_hr$Mean_HR

###Questionnaire Scoring###

##Demographics##
demos <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_demos.csv")
health <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_health.csv")
health$BMI <- (health$weight * 703)/(health$height^2)

##EPDS##

#Read EPDS Data
epds <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_epds.csv")

#Reverse Score items 3, 5-9
epds$epds_3 <- dplyr::recode(epds$epds_3, '3' = 0, '2' = 1, '1' = 2, '0' = 3)
epds$epds_5 <- dplyr::recode(epds$epds_5, '3' = 0, '2' = 1, '1' = 2, '0' = 3)
epds$epds_6 <- dplyr::recode(epds$epds_6, '3' = 0, '2' = 1, '1' = 2, '0' = 3)
epds$epds_7 <- dplyr::recode(epds$epds_7, '3' = 0, '2' = 1, '1' = 2, '0' = 3)
epds$epds_8 <- dplyr::recode(epds$epds_8, '3' = 0, '2' = 1, '1' = 2, '0' = 3)
epds$epds_9 <- dplyr::recode(epds$epds_9, '3' = 0, '2' = 1, '1' = 2, '0' = 3)

#If "Prefer Not to Answer", set to NA
epds[epds== 4] <- NA

#Score Total Depression; Mean Impute any "Prefer Not to Answer"
epds$Depression_Total <- rowMeans(epds[2:10], na.rm = TRUE) * 9 

##IAS##

#Read IAS Data
ias <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_ias.csv")

#Set "Prefer Not to Answer" to NA
ias[ias== 6] <- NA

#Calculate IAS Total Score (mean impute missing responses)
ias$IAS_Total <- rowMeans(ias[3:23], na.rm = TRUE) * 21


##MAIA-2##

#Read MAIA-2 Data
maia <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_maia2.csv")

#MAIA-2
#1. Noticing: Awareness of uncomfortable, comfortable, and neutral body sensations
#Q1______ + Q2______ + Q3______ + Q4______ / 4 = ___________
#2. Not-Distracting: Tendency not to ignore or distract oneself from sensations of pain or discomfort
#Q5(R)____ + Q6(R)____+ Q7(R)____+ Q8(R)____+Q9(R)____+Q10(R) / 6 = ___________
#3. Not-Worrying: Tendency not to worry or experience emotional distress with sensations of pain or discomfort
#Q11(R)______ + Q12(R)______ + Q13______ + Q14______ + Q15 (R) / 5 = ___________
#4. Attention Regulation: Ability to sustain and control attention to body sensations
#Q16_____ + Q17_____ + Q18_____ + Q19_____ + Q20_____ + Q21_____ + Q22_____ / 7 = ________
#5. Emotional Awareness: Awareness of the connection between body sensations and emotional states
#Q23_____ + Q24_____ + Q25_____ + Q26_____ + Q27_____ / 5 = ___________
#6. Self-Regulation: Ability to regulate distress by attention to body sensations
#Q28_____ + Q29_____ + Q30_____ + Q31_____ / 4= ___________
#7. Body Listening: Active listening to the body for insight
#Q32_____ + Q33_____ + Q34_____ / 3= ___________
#8. Trusting: Experience of one’s body as safe and trustworthy
#Q35_____ + Q36_____ + Q37_____ / 3= ___________
#Note: (R): reverse-score (5 – x) items 5, 6, 7, 8, 9 and 10 on Not-Distracting, and items 11, 12 and 15 on NotWorrying.

#Set Prefer Not to Answer to N/A
maia[maia== 6] <- NA

#Reverse Score Items 5,6,7,8,9,10,11,12,15
maia$maia_5 <- 5 - maia$maia_5
maia$maia_6 <- 5 - maia$maia_6
maia$maia_7 <- 5 - maia$maia_7
maia$maia_8 <- 5 - maia$maia_8
maia$maia_9 <- 5 - maia$maia_9
maia$maia_10 <- 5 - maia$maia_10
maia$maia_11 <- 5 - maia$maia_11
maia$maia_12 <- 5 - maia$maia_12
maia$maia_15 <- 5 - maia$maia_15

#Score all sub-scales
maia$Noticing <- rowMeans(maia[2:5], na.rm = TRUE)
maia$Not_Distracting <- rowMeans(maia[6:11], na.rm = TRUE)
maia$Not_Worrying <- rowMeans(maia[12:16], na.rm = TRUE)
maia$Attention_Regulation <- rowMeans(maia[17:23], na.rm = TRUE)
maia$Emotional_Awareness <- rowMeans(maia[24:28], na.rm = TRUE)
maia$Self_Regulation <- rowMeans(maia[29:32], na.rm = TRUE)
maia$Body_Listening <- rowMeans(maia[33:35], na.rm = TRUE)
maia$Trusting <- rowMeans(maia[36:38], na.rm = TRUE)

##QUIC##

#Read QUIC Data
quic <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_quic.csv")

#Change"Prefer Not to Answer"  to NA
quic[quic== 3] <- NA

#QUIC Scoring:
#Parental monitoring and involvement = 1R + 3R + 4R + 5R + 6R + 7R + 9R + 10R + 14R
#Parental predictability = 2 + 8R + 11 + 12 + 15R + 16 + 17R + 31 + 32 + 33 + 34 + 35
#Parental environment = 18 + 19 + 21 + 22 + 28R + 29 + 30
#Physical environment = 13 + 20 + 26 + 27 + 36R + 37 + 38
#Safety and security = 23 + 24 + 25
#Overall = Sum of all subscales.

#Fix all 'No' responses to be 0, instead of 2
quic[quic== 2] <- 0

#Reverse Score Item
quic$quic_1r <- dplyr::recode(quic$quic_1, '0' = 1, '1' = 0)
quic$quic_3r <- dplyr::recode(quic$quic_3, '0' = 1, '1' = 0)
quic$quic_4r <- dplyr::recode(quic$quic_4, '0' = 1, '1' = 0)
quic$quic_5r <- dplyr::recode(quic$quic_5, '0' = 1, '1' = 0)
quic$quic_6r <- dplyr::recode(quic$quic_6, '0' = 1, '1' = 0)
quic$quic_7r <- dplyr::recode(quic$quic_7, '0' = 1, '1' = 0)
quic$quic_8r <- dplyr::recode(quic$quic_8, '0' = 1, '1' = 0)
quic$quic_9r <- dplyr::recode(quic$quic_9, '0' = 1, '1' = 0)
quic$quic_10r <- dplyr::recode(quic$quic_10, '0' = 1, '1' = 0)
quic$quic_14r <- dplyr::recode(quic$quic_14, '0' = 1, '1' = 0)
quic$quic_15r <- dplyr::recode(quic$quic_15, '0' = 1, '1' = 0)
quic$quic_17r <- dplyr::recode(quic$quic_17, '0' = 1, '1' = 0)
quic$quic_28r <- dplyr::recode(quic$quic_28, '0' = 1, '1' = 0)
quic$quic_36r <- dplyr::recode(quic$quic_36, '0' = 1, '1' = 0)

#Sub-scales
#Parental monitoring and involvement = 1R + 3R + 4R + 5R + 6R + 7R + 9R + 10R + 14R
quic$Parental_Monitoring <- rowMeans(quic[c("quic_1r", "quic_3r","quic_4r","quic_5r",
                                            "quic_6r", "quic_7r", "quic_9r", "quic_10r", 
                                            "quic_14r")], na.rm = TRUE) * 9

#Parental predictability = 2 + 8R + 11 + 12 + 15R + 16 + 17R + 31 + 32 + 33 + 34 + 35
quic$Parental_Predictablity <- rowMeans(quic[c("quic_2", "quic_8r","quic_11","quic_12","quic_15r",
                                               "quic_16", "quic_17r", "quic_31", "quic_32", 
                                               "quic_33", "quic_34", "quic_35")], na.rm = TRUE) * 12

#Parental environment = 18 + 19 + 21 + 22 + 28R + 29 + 30
quic$Parental_Environment <- rowMeans(quic[c("quic_18", "quic_19","quic_21","quic_22","quic_28r",
                                             "quic_29", "quic_30")], na.rm = TRUE) * 7

#Physical environment = 13 + 20 + 26 + 27 + 36R + 37 + 38
quic$Physical_Environment <- rowMeans(quic[c("quic_13", "quic_20","quic_26","quic_27","quic_36r","quic_37",
                                             "quic_38")], na.rm = TRUE) * 7

#Safety and security = 23 + 24 + 25
quic$Safety <- rowMeans(quic[c("quic_23", "quic_24","quic_25")], na.rm = TRUE) * 3


#QUIC Total
quic$QUIC_Total <- quic$Parental_Monitoring + quic$Parental_Predictablity + quic$Parental_Environment +
  quic$Physical_Environment + quic$Safety

##CTQ##

#Read CTQ Data
ctq <- read.csv("/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM/In-Person_Pilot/Data/Questionnaires/isom_ctq.csv")

#Set NAs to 0 (No mean imputation for CTQ, as to not assume Trauma)
ctq[is.na(ctq)] <- 0

#Calculate Cumulative Trauma (Sum of Trauma Intensities)
ctq$Cumulative_Trauma <- ctq$ctq_1b + ctq$ctq_2b + ctq$ctq_3b + ctq$ctq_4b + ctq$ctq_5b + ctq$ctq_6c

##All Questionnaires##

#Merge all questionnaire data
isom_questionnaires <- list(epds, ias, maia, quic, ctq, demos, health)
isom_questionnaires <- isom_questionnaires %>% reduce(full_join, by='record_id')

#Fix participant_id column
colnames(isom_questionnaires)[which(names(isom_questionnaires) == "record_id")] <- "participant_id"

#Keep only Important Columns
isom_q_measures <- isom_questionnaires[c("participant_id", "participant_age", "pregnancy_week", "participant_race", "participant_ethnicity",
                                         "participant_education", "participant_income", "household_size", "participant_ses_rung",
                                         "participant_employment", "Cumulative_Trauma", "Parental_Monitoring","Parental_Predictablity",
                                         "Parental_Environment" ,"Physical_Environment","Safety","QUIC_Total",  "Depression_Total", "BMI", "current_stress", "Antidepressant", "Birth_Control",
                                         "hours_slept", "caffeine_24hr", "caffeine_time", "soda_24hr", "food_time")]

#Combine with HBD Data
isom_measures_1 <- merge(isom_q_measures, hbd.dprime, by = "participant_id")

#Add Heartrate data
isom_measures <- merge(isom_measures_1, isom_hr, by = "participant_id")

#Create Group Variable
isom_measures$Group <- ifelse(isom_measures$pregnancy_week == 0, "Comparison", "Pregnant")

#Remove isom_p36 -- too many prefer not to answers on questionnaires/low effort during session
isom_measures <- isom_measures[!(isom_measures$participant_id == "isom_p36"),]

#Set Education as factor
isom_measures$participant_education <- as.factor(isom_measures$participant_education)


#Make Separate comparison and pregnant DFs
isom_comp <- isom_measures[(isom_measures$Group == "Comparison"),]
isom_preg <- isom_measures[(isom_measures$Group == "Pregnant"),]

#Pregnancy Week Histogram
hist(isom_preg$pregnancy_week, main = "Pregnancy Week", xlab = "Week of Pregnancy")
mean(isom_preg$pregnancy_week)
sd(isom_preg$pregnancy_week)

#Correlation Matrix
isom_vars <- isom_measures[c("participant_age", "BMI", "current_stress", "hours_slept", 
                             "A.normalized","bias", "Depression_Total", "QUIC_Total", "Cumulative_Trauma",
                              "participant_ses_rung", "meanIBI", "Confidence")]

isom_corrs <- rcorr(as.matrix(isom_vars))
corrplot(isom_corrs$r, type="lower", order="hclust", 
         p.mat = isom_corrs$P, sig.level = 0.05, insig = "pch")

summary(aov(isom_measures$participant_age ~ isom_measures$participant_race))
summary(aov(isom_measures$BMI ~ isom_measures$participant_race))
summary(aov(isom_measures$current_stress ~ isom_measures$participant_race))
summary(aov(isom_measures$hours_slept ~ isom_measures$participant_race))
summary(aov(isom_measures$A.normalized ~ isom_measures$participant_race))
summary(aov(isom_measures$Depression_Total ~ isom_measures$participant_race))
summary(aov(isom_measures$Cumulative_Trauma ~ isom_measures$participant_race))
summary(aov(isom_measures$QUIC_Total ~ isom_measures$participant_race))
summary(aov(isom_measures$bias ~ isom_measures$participant_race))


summary(aov(isom_measures$participant_age ~ isom_measures$participant_education))
summary(aov(isom_measures$BMI ~ isom_measures$participant_education))
summary(aov(isom_measures$current_stress ~ isom_measures$participant_education))
summary(aov(isom_measures$hours_slept ~ isom_measures$participant_education))
summary(aov(isom_measures$A.normalized ~ isom_measures$participant_education))
summary(aov(isom_measures$Depression_Total ~ isom_measures$participant_education))
summary(aov(isom_measures$Cumulative_Trauma ~ isom_measures$participant_education))
summary(aov(isom_measures$QUIC_Total ~ isom_measures$participant_education))
summary(aov(isom_measures$bias ~ isom_measures$participant_education))


summary(aov(isom_measures$participant_age ~ isom_measures$participant_income))
summary(aov(isom_measures$BMI ~ isom_measures$participant_income))
summary(aov(isom_measures$current_stress ~ isom_measures$participant_income))
summary(aov(isom_measures$hours_slept ~ isom_measures$participant_income))
summary(aov(isom_measures$A.normalized ~ isom_measures$participant_income))
summary(aov(isom_measures$Depression_Total ~ isom_measures$participant_income))
summary(aov(isom_measures$Cumulative_Trauma ~ isom_measures$participant_income))
summary(aov(isom_measures$QUIC_Total ~ isom_measures$participant_income))
summary(aov(isom_measures$bias ~ isom_measures$participant_income))

#Group Differences
t.test(isom_measures$participant_age ~ isom_measures$Group)
t.test(isom_measures$BMI ~ isom_measures$Group)
t.test(isom_measures$current_stress ~ isom_measures$Group)
t.test(isom_measures$hours_slept ~ isom_measures$Group)
t.test(isom_measures$participant_ses_rung ~ isom_measures$Group)
t.test(isom_measures$meanIBI ~ isom_measures$Group)
t.test(isom_measures$bias ~ isom_measures$Group)
chisq.test(isom_measures$participant_education, isom_measures$Group)
chisq.test(isom_measures$participant_income, isom_measures$Group)
chisq.test(isom_measures$participant_race, isom_measures$Group)


#Below Chance HDT Performance (Chance A' Normalized)
isom_measures$HDT_Performance <- ifelse(isom_measures$A.normalized > (pi/2), "Above Chance", "Below Chance")
table(isom_measures$HDT_Performance, isom_measures$Group)
chisq.test(isom_measures$HDT_Performance, isom_measures$Group)
hist(isom_measures$A.normalized, xlab = "Interoceptive Accuracy (A' Normalized; HBD Task)", main = "HDT Performance")


#Confidence & Bias
cor.test(isom_measures$A.normalized, isom_measures$bias)
cor.test(isom_measures$A.normalized, isom_measures$Confidence)
cor.test(isom_measures$bias, isom_measures$Confidence)

t.test(isom_measures$Confidence ~ isom_measures$Group)

#QUIC Cronbach Alpha
isom_questionnaires <- isom_questionnaires[!(isom_questionnaires$participant_id == "isom_p36"),]

quic_items <- isom_questionnaires[c("quic_1r","quic_2","quic_3r","quic_4r","quic_5r","quic_6r",
                                    "quic_7r","quic_8r","quic_9r","quic_10r","quic_11","quic_12",
                                    "quic_13","quic_14r","quic_15r","quic_16","quic_17r","quic_18",
                                    "quic_19","quic_20","quic_21","quic_22","quic_23","quic_24",
                                    "quic_25","quic_26","quic_27","quic_28r","quic_29","quic_30",
                                    "quic_31","quic_32","quic_33","quic_34","quic_35","quic_36r",
                                    "quic_37","quic_38")]

psych::alpha(quic_items) #alpha = 0.9

#EPDS Cronbach Alpha
epds_items <- isom_questionnaires[c("epds_1","epds_2","epds_3","epds_4","epds_5","epds_6",
                                    "epds_7","epds_8","epds_9")]

psych::alpha(epds_items) #alpha = 0.86

#1) Group Difference in HDT Performance
summary(lm(A.normalized ~ Group + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_measures))

HBDxGroup <- ggplot(data = isom_measures, aes(x=Group, y=A.normalized))+
  geom_jitter(size = 2, aes(color=Group))+
  #geom_smooth(method = 'lm', se=F, aes(color=Group))+
  xlab("Pregnancy Status")+
  stat_summary(fun.y=mean, geom="crossbar", aes(color= Group))+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="errorbar", aes(color= Group), width=0.33)+#, shape=23, size=2)+
  ylab("Interoceptive Accuracy (A' Normalized; HDT)")+
  theme_classic()+
  theme(text = element_text(size=20))

#2) Group Difference INT w/ CTQ or QUIC
summary(lm(A.normalized ~ Group * Cumulative_Trauma + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_measures))
summary(lm(A.normalized ~ Group * QUIC_Total + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_measures))

HDTxGroupxCTQ <- ggplot(data = isom_measures, aes(x=Cumulative_Trauma, y=A.normalized))+
  geom_point(size = 2, aes(color=Group))+
  geom_smooth(method = 'lm', se=F, aes(color=Group))+
  xlab("Cumulative Trauma (CTQ)")+
  ylab("Interoceptive Accuracy (A' normalized; HDT)")+
  theme_classic()+
  theme(text = element_text(size=20))

HDTxGroupxQUIC <- ggplot(data = isom_measures, aes(x=QUIC_Total, y=A.normalized))+
  geom_point(size = 2, aes(color=Group))+
  geom_smooth(method = 'lm', se=F, aes(color=Group))+
  xlab("QUIC Total Score")+
  ylab("Interoceptive Accuracy (A' normalized; HDT)")+
  theme_classic()+
  theme(text = element_text(size=20))

#3) HBD ~ Week * Trauma or QUIC
summary(lm(A.normalized ~ pregnancy_week * Cumulative_Trauma + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_preg))
summary(lm(A.normalized ~ pregnancy_week * QUIC_Total + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_preg))

#Testing Simple Slopes
int_mod <- lm(A.normalized ~ pregnancy_week * Cumulative_Trauma + participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_preg)
sim_slopes(model = int_mod, pred = pregnancy_week, modx = Cumulative_Trauma, modx.values = c(0, 5.42, 11.25))

#Check if this is just driven by age
cor.test(isom_preg$pregnancy_week, isom_preg$Cumulative_Trauma)

#Median Split ELA for Visualization
isom_preg$Trauma_Group <- ifelse(isom_preg$Cumulative_Trauma > median(isom_preg$Cumulative_Trauma), "High", "Low")
isom_preg$QUIC_Group <- ifelse(isom_preg$QUIC_Total > median(isom_preg$QUIC_Total), "High", "Low")

HDTxWeekxCTQ <- ggplot(data = isom_preg, aes(x=pregnancy_week, y=A.normalized))+
  geom_point(size = 2, aes(color=Trauma_Group))+
  geom_smooth(method = 'lm', se=F, aes(color=Trauma_Group))+
  xlab("Week of Pregnancy")+
  ylab("Interoceptive Accuracy (A' normalized; HDT)")+
  theme_classic()+
  theme(text = element_text(size=20))

HDTxWeekxQUIC <- ggplot(data = isom_preg, aes(x=pregnancy_week, y=A.normalized))+
  geom_point(size = 2, aes(color=QUIC_Group))+
  geom_smooth(method = 'lm', se=F, aes(color=QUIC_Group))+
  xlab("Week of Pregnancy")+
  ylab("Interoceptive Accuracy (A' normalized; HDT)")+
  theme_classic()+
  theme(text = element_text(size=20))


#4) HDT & Depression
summary(lm(Depression_Total ~ A.normalized + Group+ participant_age + participant_education + meanIBI + BMI + hours_slept, data = isom_measures))

DepressionxHDT <- ggplot(data = isom_measures, aes(x=A.normalized, y=Depression_Total))+
  geom_point(size = 2, aes(color=Group))+
  geom_smooth(method = 'lm', se=F, color="black")+
  xlab("Interoceptive Accuracy (A' normalized; HDT)")+
  ylab("Depression Total Score (EPDS)")+
  theme_classic()+
  theme(text = element_text(size=20))

#Figures

#Fig2
HBDxGroup
#ggsave(filename = "Figure_2.jpg", path = "/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM", width = 15, height = 10)

#Fig3
ggarrange(HDTxGroupxCTQ,HDTxGroupxQUIC,   
          ncol = 2,nrow = 1, common.legend = T, legend = "bottom",
          labels = c("A", "B"))
#ggsave(filename = "Figure_3.jpg", path = "/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM", width = 15, height = 10)


#Fig4
ggarrange(HDTxWeekxCTQ,HDTxWeekxQUIC,   
          ncol = 2,nrow = 1, common.legend = F, legend = "bottom",
          labels = c("A", "B"))
#ggsave(filename = "Figure_4.jpg", path = "/Users/paulsavoca/Library/CloudStorage/Box-Box/iSoM", width = 15, height = 10)

