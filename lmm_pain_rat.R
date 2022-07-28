# linear mixed-effects models of pain ratings

# the following code computes linear mixed-effects models
# to examine differences between cannabis users and cannabis
# non-users regarding pain ratings.

# ==============================
# initial setup
# 1. clear environment
# 2. libraries
# 3. import raw data


# ==============================
# 1. clear environment
rm(list = ls(all.names = TRUE))


# ---------------
# 2. libraries
# install.packages('correlation')
library(correlation)
# install.packages('dplyr')
library(dplyr)
# install.packages('ggplot2')
library(ggplot2)
# install.packages('ggpubr')
library(ggpubr)
# install.packages('lme4')
library(lme4)
# install.packages('lmerTest')
library(lmerTest)
# install.packages('patchwork')
library(patchwork)
# install.packages('performance')
library(performance)
# install.packages('see')
library(see)
# install.packages('sjPlot')
library(sjPlot)
# install.packages('tidyverse')
library(tidyverse)


# ---------------
# 3. import raw data
datum = read.csv('../data/data_scored.csv',
                 stringsAsFactors = FALSE) # wide form


# ==============================
# data wrangling/cleaning
# 1. handle univar outliers
# 2. rename cols
# 3. subset data
# 4. report mean rating across trials by group

# ==============================
# 1. handle univar outliers
# cols to clean
cols_to_clean = c('S1','S2','S3','S4','S5',
                  'S6','S7','S8','S9','S10')


# locate/adjust univariate outliers 
for (col in cols_to_clean) {
  
  # select col
  var = datum[[col]]
  
  # determine quarts
  quarters = quantile(var, na.rm=TRUE)
  
  # determine IQR, 'fence' vals (2/- 2 IQRs)
  IQR = quarters[4]-quarters[2] # 75% - 25% (IQR)
  var[var>median(var)+(2*IQR)] = median(var)+(2*IQR)
  var[var<median(var)-(2*IQR)] = median(var)-(2*IQR)
  
  # print count obs above/below fences
  print(paste0('Observations above/below fences:', col))
  print(table(datum[, col] > median(datum[, col])+(2*IQR)))
  print(table(datum[, col] < median(datum[, col])-(2*IQR)))
  
  # append col to datum
  datum[,paste0(col, '.clean')] = var
  rm(col, IQR, quarters, var)
  
}

# cleanup
rm(cols_to_clean)


# ---------------
# 2. rename cols
# rename cols for analyses, plotting
names(datum)[names(datum) == 'DQ_1'] = 'age'
names(datum)[names(datum) == 'DQ_2'] = 'sex'
names(datum)[names(datum) == 'ALCHscore'] = 'alch'
names(datum)[names(datum) == 'BAIscore'] = 'BAI'
names(datum)[names(datum) == 'BPIscore'] = 'BPI'
names(datum)[names(datum) == 'MSHQ_2_1'] = 'can.30d'
names(datum)[names(datum) == 'MSHQ_6'] = 'can.onst'
names(datum)[names(datum) == 'MSHQ_8'] = 'can.yrs'
names(datum)[names(datum) == 'MSHQ_9_1'] = 'can.rcnt.mag'
names(datum)[names(datum) == 'MSHQ_10'] = 'can.rcnt.fre'


# ---------------
# 3. compute cannabis recency col
datum$can.rcnt.prod = datum$can.rcnt.mag*datum$can.rcnt.fre


# ---------------
# 4. subset data
sdatum = datum
sdatum$BPI[is.na(sdatum$BPI)] = 0
sdatum = subset(sdatum, BPI<=4)
datum = sdatum


# ---------------
# 5. report mean rating across trials by group

# select cols
tdatum = datum[, grep('S[[:digit:]]*$', colnames(datum))]

# compute means across cols
tdatum$S.mean = rowMeans(tdatum[,1:10], na.rm = TRUE)

# combine w/ datum
datum = cbind(datum, 'S.mean'=tdatum[, 11])

# cleanup
rm(tdatum)

# group means
datum %>%
  group_by(group) %>%
  dplyr::summarize(Mean = mean(S.mean, na.rm=TRUE),
                   SD = sd(S.mean, na.rm=TRUE)) %>%
  as.data.frame


# ==============================
# data preparation- mixed models
# 1. reformat data
# 2. grand mean center IVs


# ==============================
# 1. reformat data
# reshape datum (wide-to-long)
ldatum = datum %>% gather(trial, ratings.clean,
                          S1.clean:S10.clean, 
                          factor_key=TRUE)

# reset col as numeric
ldatum$trial = as.numeric(ldatum$trial)

# find/replace vals, reset col as numeric
ldatum$sex = gsub('Male', 1, ldatum$sex)
ldatum$sex = gsub('Female', 0, ldatum$sex)
ldatum$sex = as.numeric(ldatum$sex)


# ---------------
# 2. grand mean center IVs
# determine IV cols to center
cols_to_cntr = c('age', 'alch', 'group', 'sex', 'trial', 
                 'BAI', 'can.30d', 'can.onst', 'can.yrs',
                 'can.rcnt.mag', 'can.rcnt.fre', 
                 'can.rcnt.prod')

# loop thru cols
for (col in cols_to_cntr) {
  
  # compute grand mean, subtract from obs, create new col
  ldatum[, paste0(col, '.GMC')] = ldatum[[col]] - mean(ldatum[[col]],
                                                       na.rm = TRUE)
}

# cleanup
rm(col, cols_to_cntr)


# ==============================
# models
# 1. ratings model
  # rmodel: ratings ~ trial * group
# 2. check model w/ control vars
  # rmodel.A: ratings ~ trial * group + age
  # rmodel.B: ratings ~ trial * group + sex
  # rmodel.C: ratings ~ trial * group + alcohol
# 3. check model w/ cannabis interaction
  # rmodel.CAN.1: ratings ~ trial * can.30 [users]
  # rmodel.CAN.2: ratings ~ trial * can.onst [users]
  # rmodel.CAN.3: ratings ~ trial * can.yrs [users]
  # tmodel.CAN.4: ratings ~ trial * can.prod [users]



# ==============================
# 1. ratings model

# rmodel: ratings ~ trial * group
# linear mixed-effects model w/ lmer
# IVs: time, group (can user, can non-user)
# DV: ratings
rmodel = lmer(ratings.clean ~ trial.GMC *
                 group.GMC + 
                 (1|inlabID), 
               data=ldatum)

# model performance
model_performance(rmodel)
check_model(rmodel)

# export
tab_model(rmodel, 
          file='../results/ratings/rmodel.html')

# ---------------
# 2. check model w/ control vars

# linear mixed-effects model w/ lmer
# IVs: time, group (can user, can non-user), age
# DV: ratings
# excludes outliers, includes ceiling obs
rmodel.A = lmer(ratings.clean ~ trial.GMC *
                   group.GMC +
                   age.GMC + 
                   (1|inlabID), 
                 data=ldatum)

# compare nested models via liklihood ratio test
anova(rmodel.A, rmodel)

# rmodel.B: ratings ~ trial * group + sex
# linear mixed-effects model w/ lmer
# IVs: time, group (can user, can non-user), sex
# DV: ratings
# excludes outliers, includes ceiling obs
rmodel.B = lmer(ratings.clean ~ trial.GMC * 
                   group.GMC + 
                   sex.GMC + 
                   (1|inlabID), 
                 data=ldatum)

# compare nested models via liklihood ratio test
anova(rmodel.B, rmodel)

# rmodel.C: ratings ~ trial * group + alcohol
# linear mixed-effects model w/ lmer
# IVs: time, group (can user, can non-user), alcohol
# DV: ratings
# excludes outliers, includes ceiling obs
rmodel.C = lmer(ratings.clean ~ trial.GMC * 
                   group.GMC + 
                   alch.GMC + 
                   (1|inlabID), 
                 data=ldatum)

# compare nested models via liklihood ratio test
anova(rmodel.C, rmodel)

# export
tab_model(rmodel.A, rmodel.B, rmodel.C,
          file='../results/ratings/rmodel_w_cntrls.html')


# ---------------
# 3. check model w/ cannabis interactions

# rmodel.CAN.1: ratings ~ trial * can.onst [users]
rmodel.CAN.1 = lmer(ratings.clean ~ trial.GMC *
                      can.onst.GMC +
                      (1|inlabID), 
                    data=subset(ldatum, group == '1'))

# rmodel.CAN.2: ratings ~ trial * can.yrs [users]
rmodel.CAN.2 = lmer(ratings.clean ~ trial.GMC *
                      can.yrs.GMC +
                      (1|inlabID), 
                    data=subset(ldatum, group == '1'))

# rmodel.CAN.3: ratings ~ trial * can.30 [users]
rmodel.CAN.3 = lmer(ratings.clean ~ trial.GMC *
                      can.30d.GMC + 
                      (1|inlabID), 
                    data=subset(ldatum, group == '1'))

# tmodel.CAN.4: tolerance ~ trial * can.rcnt.prod [users]
rmodel.CAN.4 = lmer(ratings.clean ~ trial.GMC *
                      can.rcnt.prod.GMC +
                      (1|inlabID), 
                    data=subset(ldatum, group == '1'))


tab_model(rmodel.CAN.1, rmodel.CAN.2, 
          rmodel.CAN.3, rmodel.CAN.4)

# export
tab_model(rmodel.CAN.1, rmodel.CAN.2, 
          rmodel.CAN.3, rmodel.CAN.4,
          file='../results/ratings/rmodel_w_canusechracs.html')


# ==============================
# plotting
# 1. setup (fix cols, helper funcs)
# 2. panel w/ bar plot, group trends, individ trends


# ==============================
# 1. setup (fix cols, helper funcs)
# fix cols
colnames(ldatum) = make.unique(names(ldatum))
ldatum$ratings = as.numeric(ldatum$ratings)
ldatum$group = as.factor(ldatum$group)
ldatum$trial = as.factor(ldatum$trial)

# function to calculate mean, standard error by trial.
# adapted from http://www.sthda.com/english/wiki/ggplot2
# -error-bars-quick-start-guide-r-software-and-data-visualization
# data : some dataframe
# varname :  col to summarize
# groupnames : col(s) to group by
data_summary = function(data, varname, groupnames){
  require(plyr)
  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      len = length(x[[col]]),
      se = sd(x[[col]],na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  summary = ddply(data, groupnames, .fun=summary_func,
                  varname)
  summary = rename(summary, c('mean'=varname))
  return(summary)
}

# use function to create summary stats
df = data_summary(ldatum, varname='ratings',
                  groupnames=c('group', 'trial'))

# bar plot
bars = ggplot(df, aes(x=trial, y=ratings, fill=group)) + 
  geom_bar(stat='identity', position='dodge',
           show.legend=FALSE) +
  geom_errorbar(aes(ymin=ratings-se, ymax=ratings+se),
                width=.2, position=position_dodge(.9)) +
  coord_cartesian(ylim=c(0, 100)) +
  scale_fill_manual(values=c('gray50','#009E73')) +
  labs(x='Trial', y='Rating (VAS)') +
  theme_classic() + 
  ggtitle('Pain Ratings x Trial') +
  theme(plot.title = element_text(hjust = 0.5))

# group trends
mlines = ggplot(ldatum, aes(x=trial, y=ratings, 
                            group=group, color=group)) +
  geom_smooth(method='lm', se=TRUE, size=.5) +
  scale_color_manual(values=c('gray50','#009E73')) +
  coord_cartesian(ylim=c(0, 100)) +
  labs(x='Trial', y='Rating (VAS)') +
  theme_classic() +
  theme(legend.position='none') + 
  ggtitle('Group Trends') +
  theme(plot.title=element_text(hjust = 0.5),
        axis.title.y=element_blank())

# individual trends
lines = ggplot(ldatum, aes(x=trial, y=ratings, 
                           group=inlabID, color=group)) +
  geom_smooth(method='lm', se=FALSE, size=.5) +
  scale_color_manual(values=c('gray50','#009E73')) +
  coord_cartesian(ylim=c(0, 100)) +
  labs(x='Trial', y='Rating (VAS)') +
  theme_classic() +
  theme(legend.position='none') + 
  ggtitle('Individual Trends') +
  theme(plot.title=element_text(hjust = 0.5),
        axis.title.y=element_blank())

# combine in panel
png(file='../results/ratings/panel_trial_x_ratings.png',
    width=225, height=75, units='mm', res=300) # init plot
ggarrange(bars, mlines, lines, ncol=3, nrow=1, labels=c('A', 'B', 'C'))
dev.off() # close plot tool
