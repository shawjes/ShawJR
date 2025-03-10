library(dplyr)
library(tidyr)
library(data.table)
library(broom)
library(tibble)
library(sjstats)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(tibble)

# Identify the best model parameterization using AIC
#Goal: Choose one model that makes all analytes 'happy'

# First identify the best fitting model using fixed effects only, no interaction terms
set.seed(1234)

m0 = by_analyte %>%
  do(fit = lm(log2(CalcConcMean_AfterImputation) ~ 1, data = .))
names(m0$fit)<-m0$Panel_Analyte

m1 = by_analyte %>%
  do(fit = lm(log2(CalcConcMean_AfterImputation) ~ T21, data = .))
names(m1$fit)<-m1$Panel_Analyte

m2 = by_analyte %>%
  do(fit = lm(log2(CalcConcMean_AfterImputation) ~ T21 + Age, data = .))
names(m2$fit)<-m2$Panel_Analyte

m3 = by_analyte %>%
  do(fit = lm(log2(CalcConcMean_AfterImputation) ~ T21 + Female, data = .))
names(m3$fit)<-m3$Panel_Analyte

m4 = by_analyte %>%
  do(fit = lm(log2(CalcConcMean_AfterImputation) ~ T21 + Age + Female, data = .))
names(m4$fit)<-m4$Panel_Analyte

models<-list()
compare.AIC<-list()
compare.maxVIF2<-list()
compare.maxVIF3<-list()
compare.maxVIF4<-list()
for ( i in 1:length(m0$fit) ) {
  compare.AIC[[i]]<-AIC( m0$fit[[i]],
                         m1$fit[[i]],
                         m2$fit[[i]],
                         m3$fit[[i]],
                         m4$fit[[i]]) %>%
    rownames_to_column() %>%
    rename(Model=rowname) %>%
    mutate(Model_Call = ifelse(Model=="m0$fit[[i]]", (m0$fit[[1]]$call %>% as.character())[2],
                               ifelse(Model=="m1$fit[[i]]", (m1$fit[[1]]$call %>% as.character())[2],
                                      ifelse(Model=="m2$fit[[i]]", (m2$fit[[1]]$call %>% as.character())[2],
                                             ifelse(Model=="m3$fit[[i]]", (m3$fit[[1]]$call %>% as.character())[2],
                                                    ifelse(Model=="m4$fit[[i]]", (m4$fit[[1]]$call %>% as.character())[2], NA)))))) %>%
    mutate(Panel_Analyte = m0$Panel_Analyte[i]) %>%
    select(Panel_Analyte, Model, Model_Call, df, AIC)
  
  compare.maxVIF2[[i]]<-max(vif(m2$fit[[i]])) %>%
    as.data.frame() %>%
    `colnames<-`("Max_VIF") %>%
    mutate(Model = "m2$fit[[i]]", Panel_Analyte = m0$Panel_Analyte[i]) %>%
    select(Panel_Analyte, Model, Max_VIF)
  compare.maxVIF3[[i]]<-max(vif(m3$fit[[i]])) %>%
    as.data.frame() %>%
    `colnames<-`("Max_VIF") %>%
    mutate(Model = "m3$fit[[i]]", Panel_Analyte = m0$Panel_Analyte[i]) %>%
    select(Panel_Analyte, Model, Max_VIF)
  compare.maxVIF4[[i]]<-max(vif(m4$fit[[i]])) %>%
    as.data.frame() %>%
    `colnames<-`("Max_VIF") %>%
    mutate(Model = "m4$fit[[i]]", Panel_Analyte = m0$Panel_Analyte[i]) %>%
    select(Panel_Analyte, Model, Max_VIF)
  
}

compare.AIC_DF<-rbindlist(compare.AIC)

compare.maxVIF[[1]]

temp2<-rbindlist(compare.maxVIF2)
temp3<-rbindlist(compare.maxVIF3)
temp4<-rbindlist(compare.maxVIF4)

compare.maxVIF_DF<-rbind(temp2, temp3, temp4)

compare.AIC_DF
compare.maxVIF_DF


