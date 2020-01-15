# Code for 
# "A Systematic Review of Animal Models and Sex as a Variable in Itch Research"
# Katherine Allen
# Jan-13-2020

####  Read in data ----
# Data available upon request
dat <- read.csv("ItchData.csv", header = TRUE)

# Total counts of animal data
total_over_time_mouse <- sum(dat$Total.Mouse)
total_over_time_rat   <- sum(dat$Total.Rat)
total_over_time_NHP   <- sum(dat$Total.Non.human.Primate)
total_over_time_dog   <- sum(dat$Total.Dog)
total_articles        <- sum(dat$Animal.Model.Articles)
total_studies         <- sum(total_over_time_mouse, total_over_time_rat,
                             total_over_time_NHP, total_over_time_dog)

#### Are researchers more likely to use mice over other animal models? ----

# H_0: proportion of mice articles is less than or equal to 0.5 
# H_A: proportion of mice articles is greater than 0.5 

# If the null is not rejected, I will compare all four of the proportions. 

# Note: my n is total studies for each animal, not total papers
# as some papers study more than one animal

mice_proportion_test <- prop.test(x           = total_over_time_mouse, 
                                  n           = total_studies,
                                  alternative = "greater")

#### Are researchers more likely to use C57 background mice or CD-1? ----

# H_0: proportion of C57 background mice is less than or equal to than proportion of CD-1 mice  
# H_A: proportion of C57 background mice is greater than the proportion of CD-1 mice  

total_C57 <- sum(dat$C57Bl6.Background)
total_CD1 <- sum(dat$CD.1.or.ICR)

# note: total_over_time_mouse is not the correct value to use here for what I am computing.
# Some articles mention multiple types of mice, 
# so n is total amount of studies on each mouse. 

total_mouse_mentions <- sum(total_C57, total_CD1, dat$BALB.c, dat$Swiss.Webster, dat$NC.Nga, dat$Other.Mouse)

# Sample sizes large enough for normal approximation of the binomial
mice_type_proportion_test <- prop.test(x           = c(total_C57, total_CD1), 
                                       n           = c(total_mouse_mentions,
                                                       total_mouse_mentions), 
                                       alternative = "greater")

#### Is the change in use of mice breeds over time significant? ----

## Check C57Bl6 Background mice over time
# To test this, use a Cochran-Armitage test for trend

n    <- dat$Total.Mouse

xCDL <- dat$C57Bl6.Background

prop.trend.test(xCDL, n)

# Test ICR mice over time

xICR <- dat$CD.1.or.ICR

prop.trend.test(xICR, n)

## Test BALB mice over time

xBALB <- dat$BALB.c

prop.trend.test(xBALB, n)


## Test Swiss Webster mice over time

xSW <- dat$Swiss.Webster

prop.trend.test(xSW, n)

## Test NC Naga mice over time

xNC <- dat$NC.Nga

prop.trend.test(xNC, n)

## Test Other mice over time

xOther <- dat$Other.Mouse

prop.trend.test(xOther, n)

####  Are researchers more likely to use Sprague-Dawley rats over other strains? ----

# H_0: proportion of Sprague_Dawley rats is less than or equal to 0.5
# H_A: proportion of Sprague_Dawley rats is greater than 0.5

total_over_time_SD_rats <- sum(dat$Sprague.Dawley)

# Sample size large enough to do Normal approximation
SD_Rat_proportion_test  <- prop.test(x           = total_over_time_SD_rats, 
                                     n           = total_over_time_rat, 
                                     alternative = "greater")

####  Are researchers more likely to use M. mulatta over M. fascicularis? ----

# Sample size is not large enough to use Normal Approximation
# do an exact binomial test (could have done an exact binomial test 
# for all past analysis where prop.test was used)

# H_0: proportion of M. mulatta articles is less than or equal to proportion of  M. fascicularis 
# H_A: proportion of Macaca Mulatta articles is greater than proportion of  M. fascicularis 

total_over_time_MM_NHP <- sum(dat$Macaca.mulatta)

total_over_time_MF_NHP <- sum(dat$Macaca.facicularis)

                             
NHP_binom_test <- binom.test(x           = c(total_over_time_MM_NHP, total_over_time_MF_NHP), 
                             n           = c(total_over_time_NHP,
                                             total_over_time_NHP), 
                             alternative = "two.sided")

####  Are researchers more likely to use client owned dogs over beagles? ---- 

# H_0: proportion of client-owned dog articles is less than or equal to proportion of Beagles
# H_A: proportion of client-owned dog articles is greater than proportion of Beagles 

total_over_time_CO_dog <- sum(dat$Client.Owned)
total_over_time_B_dog  <- sum(dat$Beagle)

dog_binom_test <- binom.test(x           = c(total_over_time_CO_dog, total_over_time_B_dog), 
                             n           = c(total_over_time_dog,
                                             total_over_time_dog), 
                             alternative = "greater")

#### Are researchers more likely to report using males over other sex categories in itch research? ----
  
# H_0: proportion of articles where males are reported as used is less than or equal to 0.5 
# H_A: proportion of articles using males is is greater than 0.5 

# As with mouse models, our n needs to be total studies referencing sex, 
# not total papers

total_over_time_male <- sum(dat$Total.Male)
total_study_mentions <- sum(total_over_time_male, dat$Total.Female, dat$Total.Male.and.Female, dat$Total.Not.Specified)

sex_proportion_test <- prop.test(x           = total_over_time_male,
                                 n           = total_study_mentions,
                                 alternative = "greater")

#### When researchers report using both sexes do they consider sex as a variable? ----
  
# Sample size is small, so I use a exact binomial test

# H_0: proportion of articles where sex is NOT considered is less than or equal to the proportion of articles where sex is
# H_A: proportion of articles where sex is NOT considered is greater than the proportion of articles where sex is  

sex_not_as_variable          <- sum(dat$Total.M.vs.F.Not.Considered)
total_articles_reporting_sex <- sum(dat$Total.Male.and.Female)

sex_binom_test <- binom.test(x           = sex_not_as_variable, 
                             n           = total_articles_reporting_sex,
                             p           = 0.50, 
                             alternative = "greater")


#### Which sex category are researchers most likely to report in the different animal models? ----

## Mouse

# Check to see if there is any statistically significant difference in number of articles between the categories 

Mouse_Male         <- sum(dat$Mouse.Male.Only)
Mouse_Female       <- sum(dat$Mouse.Female.Only)
Mouse_Both         <- sum(dat$Mouse.Male.and.Female)
Mouse_NotSpecified <- sum(dat$Mouse.Not.Specified)

# Since our sample size is large enough (more than 5 in each group), we can check
# if the categories are different using a Pearson's Chi-Squared test

chi_mouse <- chisq.test(x = c(Mouse_Male, Mouse_Female, Mouse_Both, Mouse_NotSpecified))

# There is a difference. Check if male mice are used more than 50% of the time

# H_0: proportion of articles about male mice is less than or equal to 0.5 
# H_A: proportion of articles about male mice articles is greater than 0.5 

Mouse_sex_test <- prop.test(x = Mouse_Male, 
                            n = total_over_time_mouse,
                            p = 0.50, "greater")

## Rat

# H_0: proportion of articles about male rats is less than or equal to 0.5 
# H_A: proportion of articles about male rats articles is greater than 0.5 

Rat_Male         <- sum(dat$Rat.Male.Only)
Rat_Female       <- sum(dat$Rat.Female.Only)
Rat_Both         <- sum(dat$Rat.Male.and.Female)
Rat_NotSpecified <- sum(dat$Rat.Not.Specified)

# Sample Sizes are very small here, we cannot do a chi-squared or z-test

Rat_sex_test <- binom.test(x = Rat_Male, 
                           n = total_over_time_rat,
                           p = 0.50, "greater")

## Dog

Dog_Male         <- sum(dat$Dog.Male.Only)
Dog_Female       <- sum(dat$Dog.Female.Only)
Dog_Both         <- sum(dat$Dog.Male.and.Female)
Dog_NotSpecified <- sum(dat$Dog.Not.Specified)

# Sample Sizes are very small here, we cannot do a chi-squared or z-test. 

dog_sex_test <- binom.test(x           = Dog_Both, 
                           n           = total_over_time_dog,
                           p           = 0.50,
                           alternative = "greater")

## NHP

NHP_Male         <- sum(dat$NHP.Male.Only)
NHP_Female       <- sum(dat$NHP.Female.Only)
NHP_Both         <- sum(dat$NHP.Male.and.Female)
NHP_NotSpecified <- sum(dat$NHP.Not.Specified)

# Sample Sizes are very small here, we cannot do a chi-squared or z-test. 

NHP_sex_test <- binom.test(x           = NHP_Both,
                           n           = total_over_time_NHP,
                           p           = 0.5,
                           alternative = "greater")

#### When researchers used both sexes, did they consider sex as a variable? ----

## Mouse

# H_0: proportion of mice articles where sex is not considered as a variable is less than or equal to 0.5 
# H_A: proportion of mice articles where sex is not considered as a variable is greater than 0.5 
Mouse_not_sex_as_variable    <- sum(dat$Total.Mouse.M.vs.F.Not.Considered)

Mouse_sex_test <- binom.test(x = Mouse_not_sex_as_variable, 
                             n =  Mouse_Both,
                             p = 0.50, "greater")

## Rat, we do not test. Sample size is 0

## Dog
# H_0: proportion of dog articles where sex is not considered as a variable is less than or equal to 0.5 
# H_A: proportion of dog articles where sex is not considered as a variable is greater than 0.5 

Dog_not_sex_as_variable  <- sum(dat$Total.Dog.M.vs.F.Not.Considered)

Dog_sex_as_variable_test <- binom.test(x = Dog_not_sex_as_variable, 
                                       n = Dog_Both,
                                       p = 0.50, "greater")

## NHP
# H_0: proportion of NHP articles where sex is not considered as a variable is less than or equal to 0.5 
# H_A: proportion of NHP articles where sex is not considered as a variable is greater than 0.5 

NHP_not_sex_as_variable  <- sum(dat$Total.NHP.M.vs.F.Not.Considered)

NHP_sex_as_variable_test <- binom.test(x = NHP_not_sex_as_variable, 
                                       n = NHP_Both,
                                       p = 0.50, "greater")

#### Comparing pre- and post- 2015 data ----
# Was there a change in the usage of males/female/both/not specified when comparing Pre-2015 data with Post-2015 data? 
# Are males most likely to be used overall (both pre and post 2015)?

# I include 2015 in the post-2015 category

# To test this, I use a number of tests. I will explain in the code below
# Change of the usage in sexes over all
pre_2015               <- dat[dat$Year < 2015, ]
pre_2015_male          <- sum(pre_2015$Total.Male)
pre_2015_female        <- sum(pre_2015$Total.Female)
pre_2015_both          <- sum(pre_2015$Total.Male.and.Female)
pre_2015_notspecified  <- sum(pre_2015$Total.Not.Specified)
pre_2015_totalArticles <- sum(pre_2015$Animal.Model.Articles)

### Post 2015 data
post_2015               <- dat[dat$Year >= 2015, ]
post_2015_male          <- sum(post_2015$Total.Male)
post_2015_female        <- sum(post_2015$Total.Female)
post_2015_both          <- sum(post_2015$Total.Male.and.Female)
post_2015_notspecified  <- sum(post_2015$Total.Not.Specified)
post_2015_totalArticles <- sum(post_2015$Animal.Model.Articles)

# Examine specific sex categories pre and post 2015 
# Sample sizes are long enough to use the Normal Approximation to the binomial

# Is there a change in the usage of males?

male_proportion_test <- prop.test(x           = c(pre_2015_male, post_2015_male), 
                                  n           = c(pre_2015_totalArticles,
                                                  post_2015_totalArticles), 
                                  alternative = "two.sided")

# Of females?
female_proportion_test <- prop.test(x           = c(pre_2015_female, post_2015_female), 
                                    n           = c(pre_2015_totalArticles, post_2015_totalArticles), 
                                    alternative = "two.sided")

# Of both?
both_proportion_test <- prop.test(x           = c(pre_2015_both, post_2015_both), 
                                  n           = c(pre_2015_totalArticles, post_2015_totalArticles), 
                                  alternative = "two.sided")


# Of non-specified
NS_proportion_test <- prop.test(x           = c(pre_2015_notspecified, post_2015_notspecified), 
                                n           = c(pre_2015_totalArticles, post_2015_totalArticles), 
                                alternative = "two.sided")



#### Was there a change in the consideration of sex as a variable between Pre-2015 and Post 2015? ----

# To test this, I use 
# H_0: proportion of articles written prior to 2015 where sex is considered as a variable is the same 
# as the proportion of articles written during and after 2015 where sex is considered as a variable 

# H_A: proportion of articles written prior to 2015 where sex is considered as a variable 
# is NOT the same as the proportion of articles written during and after 2015 where sex is considered as a variable 

pre_2015_count_sex_as_variable   <- sum(pre_2015$Total.M.vs.F.Considered)
post_2015_count_sex_as_variable  <- sum(post_2015$Total.M.vs.F.Considered)

sex_as_variable_2015 <- binom.test(x          = c(pre_2015_count_sex_as_variable,
                                                  post_2015_count_sex_as_variable),
                                   n          = c(pre_2015_both, post_2015_both),
                                   alternative = "two.sided")

# What about overall?
# Use a Cochran-Armitage test for trend

MvsF  <- dat$Total.M.vs.F.Considered.Variable

n_sex <- dat$Total.Male.and.Female

prop.trend.test(MvsF, n_sex)

### Is there a difference in reported itch type values? ----
Acute   <- sum(dat$Acute)
Chronic <- sum(dat$Chronic)
Both    <- sum(dat$Both)

total_itch <- sum(Acute, Chronic, Both)

chisq.test(c(Acute, Chronic, Both))

## Is acute itch featured in more than 50% of itch papers? 

itch_test <- prop.test(x           = Acute, 
                       n           = total_itch,
                       alternative = "greater")
