setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/ag2bio_sb")

library(readxl)
library(sf)
library(janitor)
library(tidyverse)
library(snakecase)
library(MatchIt)
library(gtsummary)
library(grf)
library(fastDummies)
library(cobalt)

# load data extracts from earth engine
load("./data/geedat.RData")


#### Matching analysis ####

# the following data is from the sheet "Figure 2A" of the excel file 
# "Underlying data for all figures.xlsx" from the replication material.
load("./data/allspecies.RData")

# merge with gee data
allspecies <- allspecies %>%
  mutate(site = gsub("-", "_", site)) %>%
  left_join(., geedat, by = c("site" = "Unique_Station_Name")) %>%
  # filter out 10 locations with bad gps data
  filter(!is.na(accessibility)) %>%
  separate(site, into = c("site", "camera"), sep = "_") %>%
  mutate(#cluster = gsub("nonFSC|FSC", "", site),
         certification = ifelse(grepl("nonFSC", site), 0, 1)) %>%
  rename(rai = tot)


#standardise the naming convention of columns by changing capital letters and spaces to underscores
names(allspecies)<-to_snake_case(names(allspecies))
#same for camera names
#allspecies$camera<-to_snake_case(allspecies$camera)

#DATA PROCESSING:-
#Making factors, transformations etc. :
allspecies$camera <- factor(allspecies$camera)
allspecies$certification <- factor(allspecies$certification)
allspecies$pair <- factor(allspecies$cluster)
allspecies$rai <- as.numeric(allspecies$rai)
allspecies$concession <- factor(allspecies$site)

#logtransform rai
#allspecies$log_rai <- log(allspecies$rai)
#375 cameras have rai 0 and the lowest rai above 0 is 0.00427
#add fraction of an RAI to enable log transformations
allspecies$rai_plus <- (allspecies$rai+0.01)
allspecies$log_rai_plus <- log(allspecies$rai_plus)

#shorten name of log_plus transformed rai
allspecies$log_rai<- allspecies$log_rai_plus


####covariates
allspecies$visibility <- factor(allspecies$visibility)
allspecies$slope <- factor(allspecies$slope)
allspecies$fruit_trees <- factor(allspecies$fruit_trees)
allspecies$water <- factor(allspecies$water)
allspecies$elevation <- as.numeric(allspecies$elevation)
allspecies$type_of_site <- factor(allspecies$type_of_site)
allspecies$hunting <- factor(allspecies$hunting)
allspecies$dist_roads <- as.numeric(allspecies$dist_roads)
allspecies$dist_rivers <- as.numeric(allspecies$dist_rivers)
allspecies$dist_settlements <- as.numeric(allspecies$dist_settlements)
allspecies$dist_protected_areas <- as.numeric(allspecies$dist_protected_areas)


# make a matching analysis

match_vars <- c("x_0_evi", "x_10_evi", "x_15_evi", "pop_2015","accessibility", "elevation", "slope", "fruit_trees", "water",
                "visibility", "dist_settlements", "dist_roads", "dist_rivers",
                "dist_protected_areas")



# make rowsums of evi rows
allspecies$evi_5 <- rowSums(allspecies[, c("x_0_evi", "x_1_evi", "x_2_evi", "x_3_evi", "x_4_evi")])


#replicate balance table
tbl_summary(allspecies, 
            include = all_of(c("x_0_evi", "evi_5", match_vars)),
            by = certification,
            missing = "no",
            statistic = list(all_continuous() ~ "{mean} ({sd})",
                             all_categorical() ~ "{n} ({p}%)")) %>%
  add_p() %>%
  as_gt() %>%
  gt::tab_header(
    title = "Descriptive statistics of matching variables",
    subtitle = "By type of site"
  )



md1 <- MatchIt::matchit(certification ~ pop_2015 + accessibility + elevation + fruit_trees + water + dist_settlements + dist_roads + dist_rivers + dist_protected_areas,
                 data = allspecies, 
                 #caliper = 1,
                 discard = "both",
                 method = "optimal")

md2 <- MatchIt::matchit(certification ~ pop_2015 + accessibility + elevation + fruit_trees + water + dist_settlements + dist_roads + dist_rivers + dist_protected_areas,
                        data = allspecies, 
                        #caliper = 0.5,
                        exact = c("cluster"),
                        #discard = "both",
                        method = "optimal")


summary(md1)
summary(md2)



# check balance
love.plot(md1)
love.plot(md2)


mdat1 <- MatchIt::match.data(md1)
mdat2 <- MatchIt::match.data(md2)


matchmodel1 <- lmer(log_rai ~ certification + (concession|pair), 
             data = mdat1, 
             weights = weights
             ) #full model without covariates

# include unbalanced covariates
matchmodel2 <- lmer(log_rai ~ certification + pop_2015 + accessibility + elevation + (concession|pair), 
                    data = mdat1, 
                    weights = weights
) #full model without covariates

summary(matchmodel2)


# include unbalanced covariates
matchmodel3 <- lmer(log_rai ~ certification + pop_2015 + accessibility + elevation + (1|pair) + (1|concession), 
                    data = mdat1, 
                    weights = weights
) #full model without covariates
summary(matchmodel3)

# compare different models
aa <- allFit(matchmodel1)



# try with causal forest

X <- mdat1[match_vars] %>%
  fastDummies::dummy_cols(c("visibility", "slope", "fruit_trees", "water"),
                          remove_selected_columns = TRUE, 
                          remove_first_dummy = TRUE)

X <- as.matrix(X)
Y <- mdat1$log_rai
W <- as.integer(mdat1$certification)


# fit causal forest
cf <- causal_forest(X = X,
                    Y = Y,
                    W = W,
                    num.trees = 2000,
                    clusters = as.factor(mdat1$site),
                    seed = 123)

average_treatment_effect(cf)



# doubly robust estimation
forest.W <- regression_forest(X, W, tune.parameters = "all", clusters = as.factor(mdat1$site))
W.hat <- predict(forest.W)$predictions

forest.Y <- regression_forest(X, Y, tune.parameters = "all", clusters = as.factor(mdat1$site))
Y.hat <- predict(forest.Y)$predictions

forest.Y.varimp <- variable_importance(forest.Y)

# Note: Forests may have a hard time when trained on very few variables
# (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
# in selection.
selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)

tau.forest <- causal_forest(X[, selected.vars], Y, W,
                            W.hat = W.hat, Y.hat = Y.hat,
                            tune.parameters = "all",
                            clusters = as.factor(mdat1$site)
)
average_treatment_effect(tau.forest)
