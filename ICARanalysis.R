install.packages("geostan")
install.packages("SpatialEpi")
install.packages("tidybayes")

library("sf")
library("tmap")
library("spdep")
library("rstan")
library("geostan")
library("SpatialEpi")
library("tidybayes")
library("tidyverse")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("SET YOUR WORKING DIRECTORY HERE")

london_fire_data <- read.csv("london_fire_data.csv")
gl_wards <- st_read("gl_wards.shp")

#calculate expected number of fires
london_fire_data$ExpectedNum <- round(expected(population = london_fire_data$sum_hh,
                                               cases = london_fire_data$fire_count,
                                               n.strata = 1),0)

#convert spatial adjacency to nodes and edges______________
# merge the attribute table to the shapefile
spatial.data <- merge(gl_wards, london_fire_data, by = "ward_code")
#add object id
spatial.data$OBJECTID <- 1:nrow(spatial.data)

# reordering the columns
# object id, wardcode, households, firecount, expected number, factors
spatial.data <- spatial.data[, c(9,1,4,6,7,2,3,5,8)]

# turn into a spatial object
sp.object <- as(spatial.data, "Spatial")
# turn into a matrix object
adjacencyMatrix <- shape2mat(sp.object)
# extract the components for the ICAR model
extractComponents <- prep_icar_data(adjacencyMatrix)

n <- as.numeric(extractComponents$group_size)
nod1 <- extractComponents$node1
nod2 <- extractComponents$node2
n_edges <- as.numeric(extractComponents$n_edges)

# extract independent variables and store in matrix X for stan
SelectedVariables <- london_fire_data[,c(2,3,5)]
# create factor matrix
X <- model.matrix(~ 0 + imd_score + living_env + overcrowd, data = SelectedVariables)

#create dataset to be compiled in Stan
y <- spatial.data$fire_count
x <- X
e <- spatial.data$ExpectedNum

# put all components into a list object
stan.spatial.dataset <- list(N=n, N_edges=n_edges, node1=nod1, node2=nod2, Y=y,
                             X=X, offset=e, K=ncol(x))

#______________________________________________________________________________


# run stan, obtain the posterior estimation of the parameters from the model
icar_poisson_fit = stan("original_icar_poisson_model.stan", data=stan.spatial.dataset, 
                        iter=20000, control = list(max_treedepth = 12), chains=6, verbose = FALSE)

#see estimated results for alpha, beta, sigma
# remove scientific notation
options(scipen = 999)
summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))$summary

#view spatial effects phi for each area
# show first 6 rows only instead of the full 307
head(summary(icar_poisson_fit, pars=c("phi"), probs=c(0.025, 0.975))$summary)

#more output
print(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))

# diagnostic check on the rHats - insert to data frame
diagnostic.checks <- as.data.frame(summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma", "phi", "lp__"), probs=c(0.025, 0.5, 0.975))$summary)
# create binary variable
diagnostic.checks$valid <- ifelse(diagnostic.checks$Rhat < 1.1, 1, 0)
# tabulate
table(diagnostic.checks$valid)

# diagnostic = 1, so proceed with confidence

#evaluate area-specific relative risks
# show first 6 rows only instead of the full 307
head(summary(icar_poisson_fit, pars=c("rr_mu"), probs=c(0.025, 0.975))$summary)

#extract into df
# extraction key posterior results for the generated quantities 
relativeRisk.results <- as.data.frame(summary(icar_poisson_fit, pars=c("rr_mu"), probs=c(0.025, 0.975))$summary)
# now cleaning up this table up
# insert clean row numbers to new data frame
row.names(relativeRisk.results) <- 1:nrow(relativeRisk.results)
# rearrange the columns into order
relativeRisk.results <- relativeRisk.results[, c(1,4,5,7)]
# rename the columns
colnames(relativeRisk.results)[1] <- "rr"
colnames(relativeRisk.results)[2] <- "rrlower"
colnames(relativeRisk.results)[3] <- "rrupper"
colnames(relativeRisk.results)[4] <- "rHAT"

# view clean table 
head(relativeRisk.results)

#insert columns into spatial data object
# generate risk maps
# align the results to the areas in shapefile
spatial.data$rr <- relativeRisk.results[, "rr"]
spatial.data$rrlower <- relativeRisk.results[, "rrlower"]
spatial.data$rrupper <- relativeRisk.results[, "rrupper"]

#relative risk significance
# create categories to define if an area has significant increase or decrease in risk, or nothing all 
spatial.data$Significance <- NA
spatial.data$Significance[spatial.data$rrlower<1 & spatial.data$rrupper>1] <- 0    # NOT SIGNIFICANT
spatial.data$Significance[spatial.data$rrlower==1 | spatial.data$rrupper==1] <- 0  # NOT SIGNIFICANT
spatial.data$Significance[spatial.data$rrlower>1 & spatial.data$rrupper>1] <- 1    # SIGNIFICANT INCREASE
spatial.data$Significance[spatial.data$rrlower<1 & spatial.data$rrupper<1] <- -1   # SIGNIFICANT DECREASE


#cosmetics for risk map
# show distribution for risks
summary(spatial.data$rr)
hist(spatial.data$rr)

# creating the labels
RiskCategorylist <- c(">0.0 to 0.50", "0.51 to 0.75", "0.76 to 0.99", "1.00 & <1.01",
                      "1.01 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.00", 
                      "2.01 to 2.50", "2.50 to 3.00")

# create color changes for legends
RRPalette <- c("#65bafe","#98cffe","#cbe6fe","white","#fed5d5","#fcbba1","#fc9272","#fb6a4a","#de2d26","#a50f15")

# categorising the risk values to match the labelling in RiskCategorylist object
spatial.data$RelativeRiskCat <- NA
spatial.data$RelativeRiskCat[spatial.data$rr>= 0 & spatial.data$rr <= 0.50] <- -3
spatial.data$RelativeRiskCat[spatial.data$rr> 0.50 & spatial.data$rr <= 0.75] <- -2
spatial.data$RelativeRiskCat[spatial.data$rr> 0.75 & spatial.data$rr < 1] <- -1
spatial.data$RelativeRiskCat[spatial.data$rr>= 1.00 & spatial.data$rr < 1.01] <- 0
spatial.data$RelativeRiskCat[spatial.data$rr>= 1.01 & spatial.data$rr <= 1.25] <- 1
spatial.data$RelativeRiskCat[spatial.data$rr> 1.25 & spatial.data$rr <= 1.50] <- 2
spatial.data$RelativeRiskCat[spatial.data$rr> 1.50 & spatial.data$rr <= 1.75] <- 3
spatial.data$RelativeRiskCat[spatial.data$rr> 1.75 & spatial.data$rr <= 2.00] <- 4
spatial.data$RelativeRiskCat[spatial.data$rr> 2.00 & spatial.data$rr <= 2.50] <- 5
spatial.data$RelativeRiskCat[spatial.data$rr> 2.50 & spatial.data$rr <= 10] <- 6

# check to see if legend scheme is balanced
table(spatial.data$RelativeRiskCat)

#generate maps
boroughs <- st_read("London_Borough_Excluding_MHW.shp")

# map of relative risk
tm_shape(spatial.data) + 
  tm_fill("RelativeRiskCat", style = "cat", title = "Relative Risk", palette = RRPalette, labels = RiskCategorylist) +
  tm_shape(boroughs) + tm_polygons(alpha = 0, border.alpha = 0.9, border.col = "black") +
  tm_text("NAME", size = "HECTARES", legend.size.show = FALSE) +
  tm_layout(
    frame = FALSE, 
    legend.position = c("RIGHT", "BOTTOM"),
    legend.title.size = 1.2, 
    legend.text.size = 0.9) +
  tm_compass(position = c("right", "top")) + 
  tm_scale_bar(position = c("left", "bottom"))

# map of significance regions
tm_shape(spatial.data) + 
  tm_fill("Significance", style = "cat", title = "Significance Categories", 
          palette = c("#33a6fe", "white", "#fe0000"), labels = c("Significantly low", "Not Significant", "Significantly high")) +
  tm_shape(boroughs) + tm_polygons(alpha = 0, border.alpha = 0.9, border.col = "black") +
  tm_text("NAME", size = "HECTARES", legend.size.show = FALSE) +
    tm_layout(
    frame = FALSE, 
    legend.position = c("RIGHT", "BOTTOM"),
    legend.title.size = 1.2, 
    legend.text.size = 0.9,
    outer.margins = c(0.05, 0.05, 0.05, 0.05)) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("left", "bottom"))

#extract and map exceedence probabilities
# extract the exceedence probabilities from the icar_possion_fit object
# compute the probability that an area has a relative risk ratio > 1.0
threshold <- function(x){mean(x > 1.00)}
excProbrr <- icar_poisson_fit %>% spread_draws(rr_mu[i]) %>% 
  group_by(i) %>% summarise(rr_mu=threshold(rr_mu)) %>%
  pull(rr_mu)

# insert the probability exceedance values into the spatial data frame
spatial.data$excProb <- excProbrr

#cosmetics for probability exceedance map 
# create the labels for the probabilities
ProbCategorylist <- c("<0.01", "0.01-0.09", "0.10-0.19", "0.20-0.29",
                      "0.30-0.39", "0.40-0.49","0.50-0.59", "0.60-0.69",
                      "0.70-0.79", "0.80-0.89", "0.90-0.99", "1.00")

# categorising the probabilities in bands of 10s
spatial.data$ProbCat <- NA
spatial.data$ProbCat[spatial.data$excProb>=0 & spatial.data$excProb< 0.01] <- 1
spatial.data$ProbCat[spatial.data$excProb>=0.01 & spatial.data$excProb< 0.10] <- 2
spatial.data$ProbCat[spatial.data$excProb>=0.10 & spatial.data$excProb< 0.20] <- 3
spatial.data$ProbCat[spatial.data$excProb>=0.20 & spatial.data$excProb< 0.30] <- 4
spatial.data$ProbCat[spatial.data$excProb>=0.30 & spatial.data$excProb< 0.40] <- 5
spatial.data$ProbCat[spatial.data$excProb>=0.40 & spatial.data$excProb< 0.50] <- 6
spatial.data$ProbCat[spatial.data$excProb>=0.50 & spatial.data$excProb< 0.60] <- 7
spatial.data$ProbCat[spatial.data$excProb>=0.60 & spatial.data$excProb< 0.70] <- 8
spatial.data$ProbCat[spatial.data$excProb>=0.70 & spatial.data$excProb< 0.80] <- 9
spatial.data$ProbCat[spatial.data$excProb>=0.80 & spatial.data$excProb< 0.90] <- 10
spatial.data$ProbCat[spatial.data$excProb>=0.90 & spatial.data$excProb< 1.00] <- 11
spatial.data$ProbCat[spatial.data$excProb == 1.00] <- 12

# check to see if legend scheme is balanced
table(spatial.data$ProbCat)

#generate probability map output
# map of exceedance probabilities
tm_shape(spatial.data) + 
  tm_fill("ProbCat", style = "cat", title = "Probability", palette = "GnBu", labels = ProbCategorylist) +
  tm_shape(boroughs) + tm_polygons(alpha = 0, border.alpha = 0.9, border.col = "black") +
  tm_text("NAME", size = "HECTARES", legend.size.show = FALSE) +
    tm_layout(
    frame = FALSE, 
    legend.position = c("RIGHT", "BOTTOM"),
    legend.title.size = 1.2, 
    legend.text.size = 0.9) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("left", "bottom"))
