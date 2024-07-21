# check for necassary installations and load packages

if (!require("ggplot2", character.only = TRUE)) {
  install.packages("ggplot2")
}
if (!require("readxl", character.only = TRUE)) {
  install.packages("readxl")
}
if (!require("openxlsx", character.only = TRUE)) {
  install.packages("openxlsx")
}
if (!require("rstatix", character.only = TRUE)) {
  install.packages("rstatix")
}
if (!require("tidyverse", character.only = TRUE)) {
  install.packages("tidyverse")
}
if (!require("ggpubr", character.only = TRUE)) {
  install.packages("ggpubr")
}
if (!require("car", character.only = TRUE)) {
  install.packages("car")
}
if (!require("dplyr", character.only = TRUE)) {
  install.packages("dplyr")
}
if (!require("PMCMRplus", character.only = TRUE)) {
  install.packages("PMCMRplus")
}
if (!require("DescTools", character.only = TRUE)) {
  install.packages("DescTools")
}
if (!require("hrbrthemes", character.only = TRUE)) {
  install.packages("hrbrthemes")
}
if (!require("ggcorrplot", character.only = TRUE)) {
  install.packages("ggcorrplot")
}
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("caijun/ggcorrplot2")
if (!require("psych", character.only = TRUE)) {
  install.packages("psych")
}
library(ggplot2)
library(readxl)
library(openxlsx)
library(rstatix)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(car)
library(PMCMRplus)
library(DescTools)
library(hrbrthemes)
library(ggcorrplot)
library(ggcorrplot2)
library(psych)
