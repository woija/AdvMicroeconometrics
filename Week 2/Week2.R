rm(list=ls()) ## Sletter alle variable og funktioner
cat("\014") # Sletter consol teksten

setwd("/Users/simonharmat/Dropbox/Studie/Adv microeconometrics/Week 2")

# input Stata fi

library(readstata13)
mydata <- read.dta13("CreponData.dta")