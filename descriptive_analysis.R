# Descriptive analysis
# 5 jan 2026 peer review round 2

library(dplyr)
library(jtools)
library(glmtoolbox)

# Load data
wnv <- read.csv("Data/WNV seroprevalence study_Aissam_IP DZ samp dates.csv")

# Fixing test results
wnv$res <- ifelse(wnv$PRNT90.final == "Positif", 1, 0)

# Fixing location names
wnv$province <- wnv$District
wnv$province[wnv$province == "Alger"] <- "Tizi Ouzou"
wnv$province[wnv$province == "Adrar"] <- "Timimoun"

# Add age categories
wnv2 <- wnv |>
  rename(age = Age..y.) |>
  mutate(
    # Create categories
    age_group = case_when(
      age < 10             ~ "0-9",
      age >= 10 & age < 30 ~ "10-29",
      age >= 30 & age < 55 ~ "30-54",
      age >= 55 & age < 65 ~ "55-64",
      age >= 65            ~ "65+"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("0-9","10-29","30-54","55-64","65+")
    )) 

######################
# Univariable analysis
######################

## Logistic regression

# Age group
table(wnv2$age_group)
table(wnv2$age_group, wnv2$res)

model1 <- glm(res ~ age_group, 
              data = wnv2, 
              family = "binomial")

summary(model1)

round(exp(cbind(OR = coef(model1), confint(model1))), digits = 2)

# Sex
table(wnv2$sex)
wnv2$sex[wnv2$sex == ""] <- NA
wnv2$sex <- as.factor(wnv2$sex)
wnv2 <- within(wnv2, sex <- relevel(sex, ref = "M"))

table(wnv2$sex, wnv2$res)

model1 <- glm(res ~ sex, 
              data = wnv2, 
              family = "binomial")

summary(model1)

round(exp(cbind(OR = coef(model1), confint(model1))), digits = 2)

# Province residence
table(wnv2$province)

wnv2$province <- as.factor(wnv2$province)
wnv2 <- within(wnv2, province <- relevel(province, ref = "Oran"))

table(wnv2$province, wnv2$res)

model1 <- glm(res ~ province, 
              data = wnv2, 
              family = "binomial")

summary(model1)

round(exp(cbind(OR = coef(model1), confint(model1))), digits = 2)

# Region
table(wnv2$Region)
wnv2$Region <- as.factor(wnv2$Region)
wnv2 <- within(wnv2, Region <- relevel(Region, ref = "Nord"))

table(wnv2$Region, wnv2$res)

model1 <- glm(res ~ Region, 
              data = wnv2, 
              family = "binomial")

summary(model1)

round(exp(cbind(OR = coef(model1), confint(model1))), digits = 2)

# Residence type
table(wnv2$Residency.type)
wnv2$Residency.type[wnv2$Residency.type == ""] <- NA
wnv2$Residency.type <- as.factor(wnv2$Residency.type)
wnv2 <- within(wnv2, Residency.type <- relevel(Residency.type, ref = "Urbain"))

table(wnv2$Residency.type, wnv2$res)

model1 <- glm(res ~ Residency.type, 
              data = wnv2, 
              family = "binomial")

summary(model1)

round(exp(cbind(OR = coef(model1), confint(model1))), digits = 2)

# Blood transfusion
table(wnv2$Blood.transfusion)
wnv2$Blood.transfusion[wnv2$Blood.transfusion == ""] <- NA
wnv2$Blood.transfusion <- as.factor(wnv2$Blood.transfusion)
wnv2 <- within(wnv2, Blood.transfusion <- relevel(Blood.transfusion, ref = "NON"))

table(wnv2$Blood.transfusion, wnv2$res)

model1 <- glm(res ~ Blood.transfusion, 
              data = wnv2, 
              family = "binomial")

summary(model1)

exp(cbind(OR = coef(model1), confint(model1)))

# Mosquitos near living area
table(wnv2$Mousquito.risk)

wnv2$Mousquito.risk[wnv2$Mousquito.risk == ""] <- NA
wnv2$Mousquito.risk <- as.factor(wnv2$Mousquito.risk)
wnv2 <- within(wnv2, Mousquito.risk <- relevel(Mousquito.risk, ref = "OUI")) # wrong in my data apparently acc to Aissam

table(wnv2$Mousquito.risk, wnv2$res)

model1 <- glm(res ~ Mousquito.risk, 
              data = wnv2, 
              family = "binomial")

summary(model1)

exp(cbind(OR = coef(model1), confint(model1)))

# Residence near water reservoir
table(wnv2$Residency.near.water.source)

wnv2$Residency.near.water.source[wnv2$Residency.near.water.source == ""] <- NA
wnv2$Residency.near.water.source <- as.factor(wnv2$Residency.near.water.source)
wnv2 <- within(wnv2, Residency.near.water.source <- relevel(Residency.near.water.source, ref = "NON"))

table(wnv2$Residency.near.water.source, wnv2$res)

model1 <- glm(res ~ Residency.near.water.source, 
              data = wnv2, 
              family = "binomial")

summary(model1)

exp(cbind(OR = coef(model1), confint(model1)))

# Outdoor activity
table(wnv2$Actvity.type)

wnv2$Actvity.type[wnv2$Actvity.type == ""] <- NA
wnv2$Actvity.type <- as.factor(wnv2$Actvity.type)
wnv2 <- within(wnv2, Actvity.type <- relevel(Actvity.type, ref = "Indoor"))

table(wnv2$Actvity.type, wnv2$res)

model1 <- glm(res ~ Actvity.type, 
              data = wnv2, 
              family = "binomial")

summary(model1)

exp(cbind(OR = coef(model1), confint(model1)))


## Chi-square for all variables

# Province
# Create a contingency table with the needed variables.           
dat_wn <- table(wnv2$province, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Age group
dat_wn <- table(wnv2$age_group, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Sex
dat_wn <- table(wnv2$sex, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Region
dat_wn <- table(wnv2$Region, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Residence type
dat_wn <- table(wnv2$Residency.type, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Blood transfusion
dat_wn <- table(wnv2$Blood.transfusion, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Residence in mosquito infested area
dat_wn <- table(wnv2$Mousquito.risk, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Residence near water reservoir
dat_wn <- table(wnv2$Residency.near.water.source, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

# Outdoor activity
dat_wn <- table(wnv2$Actvity.type, wnv2$res) 
print(dat_wn)
print(chisq.test(dat_wn))

#########################
# Multivariable analysis
#########################
table(wnv2$Region)
# Include all variables with a p value of at least 0.2 in the univariable 
# analysis and use a stepdown procedure to remove non-significant variables
# one at a time
# complete case analysis (remove missing)
# Do not include outdoor activity here as it has more than 60% missing
wnv_complete <- wnv2 |>
  filter(!(is.na(res)) &
         !(is.na(age_group)) &
         !(is.na(province)) & 
         !(is.na(Region))  &
         !(is.na(Residency.type))  &
         !(is.na(Mousquito.risk))  & 
         !(is.na(Residency.near.water.source))) 

model2 <- glm(res ~ age_group +
                province +
                Region +
                Residency.type +
                Mousquito.risk +
                Residency.near.water.source,
              data = wnv_complete, 
              family = "binomial")

summary(model2)

# Error: 1 coefficient not defined because of singularities
alias(model2)

# Remove province
model2 <- glm(res ~ age_group +
                Region +
                Residency.type +
                Mousquito.risk +
                Residency.near.water.source,
              data = wnv_complete, 
              family = "binomial")

summary(model2)

# Remove residency near water source
model2 <- glm(res ~ age_group +
                Region +
                Residency.type +
                Mousquito.risk,
              data = wnv_complete, 
              family = "binomial")

summary(model2)

# Remove mosquito risk
model2 <- glm(res ~ age_group +
                Region +
                Residency.type,
              data = wnv_complete, 
              family = "binomial")

summary(model2)

# Re-fit with all non-missing obs for the final model
wnv_final <- wnv2 |>
  filter(!(is.na(res)) &
           !(is.na(age_group)) &
           !(is.na(province)) & 
           !(is.na(Region))  &
           !(is.na(Residency.type))) 


model2 <- glm(res ~ age_group +
                Region +
                Residency.type,
              data = wnv_final, 
              family = "binomial")

summary(model2)


exp(cbind(OR = coef(model2), confint(model2)))

table(wnv2$province, wnv2$Region)

# Check region vs province
# Which variable gives the lowest AIC?
model3 <- glm(res ~ age_group +
                province +
                Residency.type,
              data = wnv_final, 
              family = "binomial")

summary(model3)

# Model with region has lowest AIC






