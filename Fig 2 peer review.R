# Figure 2 peer review
# 19 Aug 2025
# Updated 29 March 2026 for final revisions Nat Comms

library(dplyr)
library(tidyverse)
library(ggplot2)
library(Rsero)
library(patchwork)

# Load custom plotting functions
source("custom plotting functions simple.R")
source("custom plotting functions.R")

# Load the data
wnv <- read.csv("Data/WNV seroprevalence study_Aissam_IP DZ samp dates.csv")

# Outcomes
wnv$res <- ifelse(wnv$PRNT90.final == "Positif", TRUE, FALSE)

# Remove obs with missing age
wnv <- wnv |>
  filter(!is.na(Age..y.)) |>
  rename(age = Age..y.)

# Sample times
wnv$sample_year <- as.numeric(gsub("\\D", "", wnv$Collection.date))

# Separate by province
wnv$province <- wnv$District
wnv$province[wnv$province == "Alger"] <- "Tizi Ouzou" 
wnv$province[wnv$province == "Adrar"] <- "Timimoun" 

# Provinces for which simple model was the best fit (order matters)
prov_simple_model <- c("Tizi Ouzou", "Jijel")

#Use Rsero R package to organize the data
my_data <- list()
for (i in 1:length(prov_simple_model)){
  
  prov_name <- filter(wnv, province == prov_simple_model[i])
  wnv_survey  = SeroData(age_at_sampling = prov_name$age,
                         Y               = prov_name$res, # seropositive or seronegative
                         sampling_year   = prov_name$sample_year)
  my_data[[i]] <- wnv_survey
}

# Load model fits
tizi <- readRDS("Results/chains new Rsero updated/tizi_simple.RDS")
jijel <- readRDS("Results/chains new Rsero updated/jijel_simple.RDS")
oran <- readRDS("Results/chains new Rsero updated/oran_outbreak.RDS")
tim <- readRDS("Results/chains new Rsero updated/tim_constant.RDS")
biskra <- readRDS("Results/chains new Rsero updated/biskra_constant.RDS")

# Modified plotting function
p3 <- seroprevalence.fit.custom.simple(foi_data = my_data[[1]], model_fit = tizi, colorribbon = "#B07BAC", colorline = "#694873") 
p3 <- p3[[1]] +  ggtitle("Tizi Ouzou") 

p2 <- seroprevalence.fit.custom.simple(foi_data = my_data[[2]], model_fit = jijel, colorribbon = "#78C5FC", colorline = "#027BCE") 
p2 <- p2[[1]] +  ggtitle("Jijel") 

p1 <- seroprevalence.fit.custom(FOIfit = oran, colorribbon = "#0DAB76", colorline = "#0B5D1E") 
p1 <- p1[[1]] +  ggtitle("Oran") 

p4 <- seroprevalence.fit.custom(FOIfit = tim, colorribbon = "#FFA065", colorline = "#F85E00") 
p4 <- p4[[1]] +  ggtitle("Timimoun") 

p5 <- seroprevalence.fit.custom(FOIfit = biskra, colorribbon = "#F6D370", colorline = "#E5B42D") 
p5 <- p5[[1]] +  ggtitle("Biskra") 

# Combine plots to match order in Fig 1 (map)
p1 + p3 + p2 + p4 + p5 + plot_annotation(tag_levels = "A")

# pdf 7.08 x 5 inch landscape



