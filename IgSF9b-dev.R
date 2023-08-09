library(readxl)
library(tidyverse)

setwd("~/Gao Lab Projects/IgSF9b/Dif Expression Analysis/BrainSpan/IgSF9b_dev_dlPFC")

expression <- (read.csv("Expression.csv", header = FALSE))

ages <- read.csv("Columns.csv", header = TRUE)
ages$expression <- as.numeric(expression[1,2:36]) #converted df to numeric vector 

ages$donor_age <- as.factor(ages$donor_age)

ggplot(ages, aes(x = "donor_age", y =  "expression")) +
       geom_bar()

by_expression <- ages %>% mutate(donor_age = fct_inorder(donor_age)) %>% 
  group_by(donor_age) %>% summarise(mean_exp = mean(expression))
ggplot(by_expression) + 
  geom_col(aes(x = donor_age, y = mean_exp))
