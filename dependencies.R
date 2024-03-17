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
library(ggplot2)
library(readxl)
library(openxlsx)
library(rstatix)
library(tidyverse)
library(ggpubr)
library(dplyr)