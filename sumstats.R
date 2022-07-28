# summary stats for study 1: pain ratings/tolerance

# the following code computes summary statistics
# across variables for cannabis users and cannabis
# non-users.

# ==============================
# initial setup
# 1. libraries
# 2. directories
# 3. raw data
# ==============================

# clear environment
rm(list = ls(all.names = TRUE))

# libraries
# install.packages('dplyr') # run once
library(dplyr)

# import raw data
datum = read.csv('../data/data_scored.csv') # wide form

# ==============================
# summary statistics
# 1. name cols to summ
# 2. compute summ statistics
# ==============================

# rename cols for analyses
names(datum)[names(datum) == 'DQ_1'] = 'age'
names(datum)[names(datum) == 'DQ_2'] = 'sex'
names(datum)[names(datum) == 'DQ_3'] = 'race'
names(datum)[names(datum) == 'ALCHscore'] = 'alch'
names(datum)[names(datum) == 'BPIscore'] = 'BPI'
names(datum)[names(datum) == 'MSHQ_2_1'] = 'can.30'
names(datum)[names(datum) == 'MSHQ_6'] = 'can.onst'
names(datum)[names(datum) == 'MSHQ_8'] = 'can.yrs'
names(datum)[names(datum) == 'MSHQ_9_1'] = 'can.rcnt.mag'
names(datum)[names(datum) == 'MSHQ_10'] = 'can.rcnt.fre'

# 3. compute cannabis recency col
datum$can.rcnt.prod = datum$can.rcnt.mag*datum$can.rcnt.fre

# subset data
datum$BPI[is.na(datum$BPI)] = 0
datum = subset(datum, BPI<=4)

# cols to compute summaries
summ_cols = c('age', 'height', 'weight','BMI', 
              'BDIscore', 'BAIscore', 'PHQsomscore',
              'PQBtotscore', 'can.30', 'can.onst', 
              'can.yrs', 'C48H', 'alch', 'can.rcnt.prod')

# compute summaries
summs = datum %>%
  group_by(group) %>%
  summarise_at(vars(summ_cols), funs(mean(., na.rm=TRUE),
                                     sd(., na.rm=TRUE)))
# export
write.csv(summs, '../results/sumstats.csv')

