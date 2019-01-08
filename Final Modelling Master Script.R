#########################################################################################################
##################################################### Settings ##########################################
#########################################################################################################
setwd("~/UNIL/Master/2th semester/Master Project/Data")

library(xlsx)
library(PresenceAbsence)
library(ecospat)
library(AUC)
library(DescTools)
library(devtools)
library("HMSC")
library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(randomForest)
library(gtable)
library(corrplot)

Plants_df <- read.csv("Total_sp_env_dataframe_good.csv", sep = ";")
sp_richness <- vector(mode = "numeric", length = nrow(Plants_df))
sp_richness <- apply(Plants_df[,10:109], MARGIN = 1, sum)
Plants_df <- cbind(Plants_df, sp_richness)
# list of species names
sp_names <- colnames(Plants_df)[10:109]

#########################################################################################################
##################################################### scaling the predictors ############################
#########################################################################################################
Plants_df_scaled <- Plants_df
Plants_df_scaled[,4:8] <- sapply(Plants_df_scaled[,4:8], function(x) scale(x))

#########################################################################################################
################### Adding random var distributions for species used as predictors in GLMS ##############
#########################################################################################################
sp_names[59]
length(Plants_df_scaled$Lotus_corniculatus_aggr[Plants_df_scaled$Lotus_corniculatus_aggr == "1"])
length(Plants_df_scaled$Lotus_corniculatus_aggr[Plants_df_scaled$Lotus_corniculatus_aggr == "0"])

sp_names[58]
length(Plants_df_scaled$Lolium_perenne[Plants_df_scaled$Lolium_perenne == "1"])
length(Plants_df_scaled$Lolium_perenne[Plants_df_scaled$Lolium_perenne == "0"])

sp_names[49]
length(Plants_df_scaled$Holcus_lanatus[Plants_df_scaled$Holcus_lanatus == "1"])
length(Plants_df_scaled$Holcus_lanatus[Plants_df_scaled$Holcus_lanatus == "0"])

sp_names[39]
length(Plants_df_scaled$Festuca_rubra_aggr[Plants_df_scaled$Festuca_rubra_aggr == "1"])
length(Plants_df_scaled$Festuca_rubra_aggr[Plants_df_scaled$Festuca_rubra_aggr == "0"])

sp_names[94]
length(Plants_df_scaled$Trifolium_pratense_sl[Plants_df_scaled$Trifolium_pratense_sl == "1"])
length(Plants_df_scaled$Trifolium_pratense_sl[Plants_df_scaled$Trifolium_pratense_sl == "0"])

sp_names[7]
length(Plants_df_scaled$Anthoxanthum_odoratum_aggr[Plants_df_scaled$Anthoxanthum_odoratum_aggr == "1"])
length(Plants_df_scaled$Anthoxanthum_odoratum_aggr[Plants_df_scaled$Anthoxanthum_odoratum_aggr == "0"])

Lotus_corniculatus_aggr_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(381/911,530/911))
Lolium_perenne_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(113/911,798/911))
Holcus_lanatus_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(122/911,789/911))
Festuca_rubra_agrr_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(434/911,477/911))
Trifolium_pratense_sl_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(421/911,490/911))
Anthoxanthum_odoratum_aggr_random <- sample(c(1,0), size = 911, replace = TRUE, prob = c(384/911,527/911))

Plants_df_scaled <- cbind(Plants_df_scaled,Lotus_corniculatus_aggr_random, Lolium_perenne_random,
                          Holcus_lanatus_random, Festuca_rubra_agrr_random,
                          Trifolium_pratense_sl_random, Anthoxanthum_odoratum_aggr_random)

#########################################################################################################
##################################################### 4 D Matrix ##########################################
#########################################################################################################

# Creating the results matrix
dimensions <- c(911,100,10,10)
dim_names <- c("Eval Plots", "Species", "Models", "Trials") # Dimensions content
forDMatrix <- array(data = NA, dim = dimensions)
dimnames(forDMatrix)[[1]] <- as.list(Plants_df$Plot)
dimnames(forDMatrix)[[2]] <- as.list(sp_names)
dimnames(forDMatrix)[[3]] <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4", "Random Forest corr", "Random Forest prev", "Random Forest sp richness", "GLM5")
dimnames(forDMatrix)[[4]] <- as.list(seq(1,10,1))

#########################################################################################################
##################################################### Stratified Random sampling  #######################
#########################################################################################################


# Creating a plot score to be sure that it's not always the same plots
# taken into the evaluating dataset during the random sampling
plots_scores <- as.data.frame(cbind(Plants_df$Plot, Plants_df$DEM_vec ,rep(0,911)))
colnames(plots_scores) <- c("Plot", "DEM", "score")
plots_scores <- plots_scores[order(plots_scores$DEM),]

min(Plants_df$DEM_vec,na.rm = T) #418.8
max(Plants_df$DEM_vec,na.rm = T) #3101.4
#The aim of the random sampling : take 80% of the plot for training and 20% for evaluating each 200m
elev <- seq(400,3100, by = 200)

#Creating a categorical list with all plots within 200m in one slot of a list
cate_elev <- list()
for (elevation in 1:length(elev)) {
  cate_elev[[elevation]] <- plots_scores$Plot[plots_scores$DEM >= as.numeric(elev[elevation]) &
                                                plots_scores$DEM <= as.numeric(elev[elevation]+200)]
}

#Merging the two last slots of the list (to low number of plots in those slots)
cate_elev[[13]] <- c(cate_elev[[13]], cate_elev[[14]])
cate_elev <- cate_elev[1:13]

#Creating empty lists to fill up with the results !
training_plots <- list()
evaluating_plots <- vector()
species_fitting <- list(species = vector(length = 100))
trials_fitting <- list(trial1 = species_fitting, trial2 = species_fitting, trial3 = species_fitting, trial4 = species_fitting, trial5 = species_fitting, trial6 = species_fitting, trial7 = species_fitting, trial8 = species_fitting, trial9 = species_fitting, trial10 = species_fitting)
models_fitting <- list(GLM1 = trials_fitting, GLM2 = trials_fitting, HMSC_uncond = trials_fitting, HMSC_cond = trials_fitting, GLM3 = trials_fitting, GLM4 = trials_fitting, Random_Forest_corr = trials_fitting, Random_Forest_prev = trials_fitting, Random_Forest_SR = trials_fitting, GLM5 = trials_fitting)

for (trials in 1:length(forDMatrix[1,1,1,])) {
  
  # Because at trial = 1, all plots have the same score
  if(trials == 1)
    # Looping through all categories of elevation
    for (categories in 1:length(cate_elev)) {
      plots <- unlist(cate_elev[categories])
      training_plots[[categories]] <- sample(x = plots, size = as.numeric(round(0.8*length(plots))))
    }
  
  # if all plots does not have the same score, we need first to take those with lower score
  if(trials > 1)
    # Looping through all categories of elevation
    for (categories in 1:length(cate_elev)) {
      plots <- unlist(cate_elev[categories])
      # If les than 80% of plots have the lowest score, pich them up and then complete the 80%
      if((length(plots_scores$Plot[plots_scores$Plot %in% plots & plots_scores$score == min(plots_scores$score[plots_scores$Plot %in% plots])])/length(plots))<0.8)
        # First taking the plots with lowest score
        training_plots[[categories]] <- plots_scores$Plot[plots_scores$Plot %in% plots & plots_scores$score == min(plots_scores$score[plots_scores$Plot %in% plots])]
        # Then filling up the 80% needed with random sampling among the remaining plots
        training_plots[[categories]] <- append(training_plots[[categories]],
                                                sample(x = plots[!plots %in% training_plots[[categories]]],
                                                      size = as.numeric((0.8-(length(training_plots[[categories]])
                                                                            / length(plots)))*length(plots))))
      # If more than 80% of the plots have the lowest score, pick up training plots inside them
      if((length(plots_scores$Plot[plots_scores$Plot %in% plots & plots_scores$score == min(plots_scores$score[plots_scores$Plot %in% plots])])/length(plots))>=0.8)
        training_plots[[categories]] <- sample(x = plots_scores$Plot[plots_scores$Plot %in% plots & plots_scores$score == min(plots_scores$score[plots_scores$Plot %in% plots])], size = as.numeric(round(0.8*length(plots))))
    }
        
  # Making score +1 for the plots selected in the training dataset
  plots_scores$score[plots_scores$Plot %in% unlist(training_plots)] <- 
    as.numeric(plots_scores$score[plots_scores$Plot %in% unlist(training_plots)]) + 1
  # Filling the evaluating dataset with the plots not taking in the evaluating dataset !
  evaluating_plots <- plots_scores$Plot[!plots_scores$Plot %in% unlist(training_plots)]
  
  training_Tmp <- as.data.frame(Plants_df_scaled[Plants_df_scaled$Plot %in% unlist(training_plots),])
  evaluating_Tmp <- as.data.frame(Plants_df_scaled[Plants_df_scaled$Plot %in% evaluating_plots,])
  
  #########################################################################################################
  ##################################################### Modelling GLM1 ########################################
  #########################################################################################################
  
  formula_GLM1 <- vector(length = 100)
  for(species in 1:length(formula_GLM1)){
    formula_GLM1[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + Lotus_corniculatus_aggr_random + Lolium_perenne_random + Holcus_lanatus_random", sep = "")
  }
  
  fit_GLM1 <- vector(mode = "list", length = 100)
  prediction_GLM1 <- vector(mode = "list", length = 100)
  for(species in 1:length(fit_GLM1)){
    fit_GLM1[[species]] <- glm(formula_GLM1[species], family = "binomial", data = training_Tmp)
    prediction_GLM1[[species]] <- predict(fit_GLM1[[species]], newdata = evaluating_Tmp, type = "response")
  }
  
  # Save model for Rsquared
  for (species in 1:100) {
    models_fitting$GLM1[[trials]][[species]] <- fit_GLM1[[species]]
  }
 
  #########################################################################################################
  ################################################ Fill 4D matrix GLM1 ####################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,1,trials] <- prediction_GLM1[[species]][which(names(prediction_GLM1[[species]]) == plots)]
      }
  }
  
  #########################################################################################################
  ##################################################### Modelling GLM2 ########################################
  #########################################################################################################
  
  formula_GLM2 <- vector(length = 100)
  for(species in 1:length(formula_GLM2)){
    formula_GLM2[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + Festuca_rubra_agrr_random + Trifolium_pratense_sl_random + Anthoxanthum_odoratum_aggr_random", sep = "")
  }
  
  fit_GLM2 <- vector(mode = "list", length = 100)
  prediction_GLM2 <- vector(mode = "list", length = 100)
  for(species in 1:length(fit_GLM2)){
    fit_GLM2[[species]] <- glm(formula_GLM2[species], family = "binomial", data = training_Tmp)
    prediction_GLM2[[species]] <- predict(fit_GLM2[[species]], newdata = evaluating_Tmp, type = "response")
  }
  
  # Save model for Rsquared
  for (species in 1:100) {
    models_fitting$GLM2[[trials]][[species]] <- fit_GLM2[[species]]
  }
  
  #########################################################################################################
  ############################################### Fill 4D matrix GLM 2 ####################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
    for (species in 1:100) {
      forDMatrix[plots,species,2,trials] <- prediction_GLM2[[species]][which(names(prediction_GLM2[[species]]) == plots)]
    }
  }

  #########################################################################################################
  ################################### HMSC matrices (training and evaluating) preparation #################
  #########################################################################################################
  
  # TRAINING MATRICES
  
  # Training plants matrix
  train_plant_Tmp <- training_Tmp[,10:109]
  train_plant_Tmp <- as.matrix(train_plant_Tmp)
  
  # Training ecological variables matrix
  train_var_Tmp <- training_Tmp[,4:8]
  # Adding the SECOND TERM ORDER
  second_term_mat_train <- matrix(ncol = 5, nrow = nrow(training_Tmp))
  second_term_mat_train <- as.data.frame(second_term_mat_train)
  second_term_mat_train <- sapply(train_var_Tmp, function(x) x^2)
  colnames(second_term_mat_train) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
  train_var_Tmp <- cbind(train_var_Tmp, second_term_mat_train)
  train_var_Tmp <- as.matrix(train_var_Tmp)
  
  # EVALUATING MATRICES
  
  # Evaluating plants matrix
  eval_plant_Tmp <- evaluating_Tmp[,10:109]
  eval_plant_Tmp <- as.matrix(eval_plant_Tmp)
  
  # Evaluating ecological variables matrix
  eval_var_Tmp <- evaluating_Tmp[,4:8]
  # Adding the SECOND TERM ORDER
  second_term_mat_eval <- matrix(ncol = 5, nrow = nrow(evaluating_Tmp))
  second_term_mat_eval <- as.data.frame(second_term_mat_eval)
  second_term_mat_eval <- sapply(eval_var_Tmp, function(x) x^2)
  colnames(second_term_mat_eval) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
  eval_var_Tmp <- cbind(eval_var_Tmp, second_term_mat_eval)
  eval_var_Tmp <- as.matrix(eval_var_Tmp)
  
  #########################################################################################################
  ############################################ HMSC unconditional #########################################
  #########################################################################################################
  
  # Formatting data
  Form_uncond_HMSC_train <- as.HMSCdata(Y = train_plant_Tmp, X = train_var_Tmp, scaleX = FALSE, interceptX = TRUE)
  Form_uncond_HMSC_eval <- as.HMSCdata(Y = eval_plant_Tmp, X = eval_var_Tmp, scaleX = FALSE, interceptX = TRUE)
  
  # Fitting the model
  model_uncond_HMSC <- hmsc(Form_uncond_HMSC_train, family = "probit", niter = 10000, nburn = 5000, thin = 50)
 
  # Predictions
  pred_uncond_HMSC <- predict(model_uncond_HMSC, newdata = Form_uncond_HMSC_eval, type = "response")
  
  # Save model fitting for Rsquared
  models_fitting$HMSC_uncond[[trials]]$species <- model_uncond_HMSC
  
  #########################################################################################################
  ################################################# Fill 4D matrix HMSC uncond ############################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,3,trials] <- pred_uncond_HMSC[which(rownames(pred_uncond_HMSC) == plots), species]
      }
  }
  
  #########################################################################################################
  ##################################### HMSC conditional creating RANDOM VARIABLES ########################
  #########################################################################################################
  
  # Creating a matrix for random factors
  
  Random_factor <- matrix(ncol = 2, nrow = 911)
  rownames(Random_factor) <- Plants_df_scaled$Plot
  colnames(Random_factor) <- c("Random_1", "Random_elev")
  
  # Filling the Random_factor$Random_1 column only with 1 for each plots (as they are all in the same
  # area and all sampled the same way)
  
  Random_factor[,1] <- rep(1,911)
  
  #Filling the Random_factor$Random_elev according to the elvation class (Low or high)
  for(rows in 1:nrow(Plants_df_scaled)){
    if(Plants_df_scaled$DEM_vec[rows] < median(Plants_df_scaled$DEM_vec))
      Random_factor[rows,2] <- 2
    else
      Random_factor[rows,2] <- 3
  }
  
  # Changing it into a dataframe and setting columns as factor
  Random_factor <- as.data.frame(Random_factor)
  Random_factor$Random_1 <- as.factor(Random_factor$Random_1)
  Random_factor$Random_elev <- as.factor(Random_factor$Random_elev)
  
  # Creating Random factor training and evaluating dataframes
  Random_factor_train <- Random_factor[which(rownames(Random_factor) %in% unlist(training_plots)),]
  Random_factor_eval <- Random_factor[which(rownames(Random_factor) %in% evaluating_plots),]
  
  #########################################################################################################
  ##################################### HMSC conditional building the output matrix #######################
  #########################################################################################################
  
  # Building the matrix
  dimensions_predict_HMSC_cond <- c(nrow(evaluating_Tmp), 100, 100)
  dinames_predict_HMSC_cond <- list("Plots", "Species", "Nsample")
  predictions_HMSC_cond <- array(data = NA, dim = dimensions_predict_HMSC_cond)
  
  # Setting dimensions names
  dimnames(predictions_HMSC_cond)[[1]] <- as.list(evaluating_Tmp$Plot)
  dimnames(predictions_HMSC_cond)[[2]] <- sp_names
  dimnames(predictions_HMSC_cond)[[3]] <- seq(1,100,1)
  
  #########################################################################################################
  ############################################ HMSC conditional All species ###############################
  #########################################################################################################
  
  
  # Formatting data
  Form_cond_HMSC_train <- as.HMSCdata(Y = train_plant_Tmp, X = train_var_Tmp, Random = Random_factor_train, interceptX = TRUE, scaleX = FALSE)
  Form_cond_HMSC_eval <- as.HMSCdata(Y = eval_plant_Tmp, X = eval_var_Tmp, Random = Random_factor_eval, interceptX = TRUE, scaleX = FALSE)
  
  # Fitting the model
  model_cond_HMSC <- hmsc(Form_cond_HMSC_train, family = "probit", niter = 10000, nburn = 5000, thin = 50)
  
  # Building predictions
  for(species in 1:100){
    Eval_data_Tmp <- eval_plant_Tmp
    Eval_data_Tmp[,species] <- NA
    Form_cond_HMSC_eval <- as.HMSCdata(Y = Eval_data_Tmp, X = eval_var_Tmp, Random = Random_factor_eval, interceptX = TRUE, scaleX = FALSE)
    predictions_Tmp <- predict(model_cond_HMSC, newdata = Form_cond_HMSC_eval, type = "response", conditional =
                                 colnames(Form_cond_HMSC_eval$Y), nsample = 100)[,species,]
    for (nsample in 1:100) {
      predictions_HMSC_cond[,species,nsample] <- predictions_Tmp[,nsample]
    }
  }
  
  # Save model fitting for Rsquared
  models_fitting$HMSC_cond[[trials]]$species <- model_cond_HMSC
  
  #########################################################################################################
  ############################################ Filling 4D matrix HMSC cond ################################
  #########################################################################################################
  
  # Creating a matrix for the mean for the 100 nsamples
  mean_predictions_HMSC_cond <- matrix(ncol = 100, nrow = nrow(evaluating_Tmp))
  colnames(mean_predictions_HMSC_cond) <- sp_names
  rownames(mean_predictions_HMSC_cond) <- evaluating_Tmp$Plot
  
  # Filling this matrix
  for (species in 1:100) {
    mean_predictions_HMSC_cond[,species] <- apply(predictions_HMSC_cond[,species,], 1, function(x) mean(x))
  }
  
  # Filling 4D matrix
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,4,trials] <- mean_predictions_HMSC_cond[which(rownames(mean_predictions_HMSC_cond) == plots), species]
      }
  }
  
  #########################################################################################################
  ###################################### GLM3 correlation sp as predictors ################################
  #########################################################################################################
  
  formula_GLM3 <- vector(length = 100)
  for(species in 1:length(formula_GLM3)){
    if(species == 32)
      formula_GLM3[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[73], ") + as.factor(", sp_names[38], ")", sep = "")
    else if (species == 73)
      formula_GLM3[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[32], ") + as.factor(", sp_names[38], ")", sep = "")
    else if (species == 38)
      formula_GLM3[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[32], ") + as.factor(", sp_names[73], ")", sep = "")
    else
      formula_GLM3[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[32], ") + as.factor(", sp_names[73], ") + as.factor(", sp_names[38], ")", sep = "")
  }
  
  fit_GLM3 <- vector(mode = "list", length = 100)
  prediction_GLM3 <- vector(mode = "list", length = 100)
  for(j in 1:length(fit_GLM3)){
    fit_GLM3[[j]] <- glm(formula_GLM3[j], family = "binomial", data = training_Tmp)
    prediction_GLM3[[j]] <- predict(fit_GLM3[[j]], newdata = evaluating_Tmp, type = "response")
  }
  
  # Save model for Rsquared
  for (species in 1:100) {
    models_fitting$GLM3[[trials]][[species]] <- fit_GLM3[[species]]
  }
  
  #########################################################################################################
  ########################################## Filling 4 D matrix with GLM3 #################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,5,trials] <- prediction_GLM3[[species]][which(names(prediction_GLM3[[species]]) == plots)]
      }
  }
  
  #########################################################################################################
  ###################################### GLM4 most prevalent species as predictors ########################
  #########################################################################################################
  
  formula_GLM4 <- vector(length = 100)
  for(species in 1:length(formula_GLM4)){
    if(species == 39)
      formula_GLM4[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[94], ") + as.factor(", sp_names[7], ")", sep = "")
    else if(species == 94)
      formula_GLM4[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[39], ") + as.factor(", sp_names[7], ")", sep = "")
    else if(species == 7)
      formula_GLM4[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[39], ") + as.factor(", sp_names[94], ")", sep = "")
    else
      formula_GLM4[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + as.factor(", sp_names[39], ") + as.factor(", sp_names[94], ") + as.factor(", sp_names[7], ")", sep = "")
  }
  
  fit_GLM4 <- vector(mode = "list", length = 100)
  prediction_GLM4 <- vector(mode = "list", length = 100)
  for(j in 1:length(fit_GLM4)){
    fit_GLM4[[j]] <- glm(formula_GLM4[j], family = "binomial", data = training_Tmp)
    prediction_GLM4[[j]] <- predict(fit_GLM4[[j]], newdata = evaluating_Tmp, type = "response")
  }
  
  # Save model for Rsquared
  for (species in 1:100) {
    models_fitting$GLM4[[trials]][[species]] <- fit_GLM4[[species]]
  }
  
  #########################################################################################################
  ########################################## Filling 4 D matrix with GLM4 #################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,6,trials] <- prediction_GLM4[[species]][which(names(prediction_GLM4[[species]]) == plots)]
      }
  }
  
  #########################################################################################################
  ########################################### RANDOM FOREST CORR ##########################################
  #########################################################################################################
  
  
  Species_vec <- vector(length = nrow(training_Tmp))
  
  Predictor_df_train <- training_Tmp[,c(4:8,48,103,16)]
  Predictor_df_train <- as.data.frame(Predictor_df_train)
  colnames(Predictor_df_train) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope",
                                    "Festuca_rubra_aggr", "Trifolium_pratense_sl",
                                    "Anthoxanthum_odoratum_aggr")
  
  Predictor_df_eval <- evaluating_Tmp[,c(4:8,48,103,16)]
  Predictor_df_eval <- as.data.frame(Predictor_df_eval)
  colnames(Predictor_df_eval) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope",
                                   "Festuca_rubra_aggr", "Trifolium_pratense_sl",
                                   "Anthoxanthum_odoratum_aggr")
  
  
  fit_RF_corr <- vector(mode = "list", length = 100)
  prediction_RF_corr <- vector(mode = "list", length = 100)
  
  for(species in 1:length(fit_RF_corr)){
    
    Species_vec <- training_Tmp[,colnames(training_Tmp) == sp_names[species]]
    
    if(species == 39){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[39]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[39]]}
    else if(species == 94){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[94]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[94]]}
    else if(species == 7){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[7]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[7]]}
    else{
      Predictor_df_train_loop <- Predictor_df_train
      Predictor_df_eval_loop <- Predictor_df_eval}
    
    dta_Tmp <- cbind(Predictor_df_train_loop, Species_vec)
    
    fit_RF_corr[[species]] <- randomForest(as.factor(Species_vec)~., data = dta_Tmp, ntree= 10000)
    prediction_RF_corr[[species]] <- predict(fit_RF_corr[[species]], newdata = Predictor_df_eval_loop, type = "prob")
  }
  
  for(species in 1:100){
    models_fitting$Random_Forest_corr[[trials]][[species]] <- fit_RF_corr[[species]]
  }
  
  #########################################################################################################
  ######################################## Filling 4 D matrix with RF corr ################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,7,trials] <- prediction_RF_corr[[species]][which(rownames(prediction_RF_corr[[species]]) == plots),2]
      }
  }
  
  #########################################################################################################
  ############################################ RANDOM FOREST PREV #########################################
  #########################################################################################################
  
  
  Species_vec <- vector(length = nrow(training_Tmp))
  
  Predictor_df_train <- training_Tmp[,c(4:8,67,68,58)]
  Predictor_df_train <- as.data.frame(Predictor_df_train)
  colnames(Predictor_df_train) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope",
                                    "Lolium_perenne", "Lotus_corniculatus_aggr",
                                    "Holcus_lanatus")
  
  Predictor_df_eval <- evaluating_Tmp[,c(4:8,67,68,58)]
  Predictor_df_eval <- as.data.frame(Predictor_df_eval)
  colnames(Predictor_df_eval) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope",
                                   "Lolium_perenne", "Lotus_corniculatus_aggr",
                                   "Holcus_lanatus")
  
  
  fit_RF_prev <- vector(mode = "list", length = 100)
  prediction_RF_prev <- vector(mode = "list", length = 100)
  
  for(species in 1:length(fit_RF_prev)){
    
    Species_vec <- training_Tmp[,colnames(training_Tmp) == sp_names[species]]
    
    if(species == 58){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[58]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[58]]}
    else if(species == 59){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[59]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[59]]}
    else if(species == 49){
      Predictor_df_train_loop <- Predictor_df_train[,colnames(Predictor_df_train) != sp_names[49]]
      Predictor_df_eval_loop <- Predictor_df_eval[,colnames(Predictor_df_eval) != sp_names[49]]}
    else{
      Predictor_df_train_loop <- Predictor_df_train
      Predictor_df_eval_loop <- Predictor_df_eval}
    
    dta_Tmp <- cbind(Predictor_df_train_loop, Species_vec)
    
    fit_RF_prev[[species]] <- randomForest(as.factor(Species_vec)~., data = dta_Tmp, ntree= 10000)
    prediction_RF_prev[[species]] <- predict(fit_RF_prev[[species]], newdata = Predictor_df_eval_loop, type = "prob")
  }
  
  for(species in 1:100){
    models_fitting$Random_Forest_prev[[trials]][[species]] <- fit_RF_prev[[species]]
  }
  
  #########################################################################################################
  ###################################### Filling 4 D matrix with RF prev ##################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,8,trials] <- prediction_RF_prev[[species]][which(rownames(prediction_RF_prev[[species]]) == plots),2]
      }
  }
  
  #########################################################################################################
  ############################################# RANDOM FOREST SR ##########################################
  #########################################################################################################
  
  
  Species_vec <- vector(length = nrow(training_Tmp))
  
  Predictor_df_train <- training_Tmp[,c(4:8,110)]
  Predictor_df_train <- as.data.frame(Predictor_df_train)
  colnames(Predictor_df_train) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope", "SR")
  
  Predictor_df_eval <- evaluating_Tmp[,c(4:8,110)]
  Predictor_df_eval <- as.data.frame(Predictor_df_eval)
  colnames(Predictor_df_eval) <- c("GS_S_Rad", "GS_p", "GS_t", "pH", "slope", "SR")
  
  
  fit_RF_SR <- vector(mode = "list", length = 100)
  prediction_RF_SR <- vector(mode = "list", length = 100)
  
  for(species in 1:length(fit_RF_corr)){
    
    Species_vec <- training_Tmp[,colnames(training_Tmp) == sp_names[species]]
    dta_Tmp <- cbind(Predictor_df_train, Species_vec)
    
    fit_RF_SR[[species]] <- randomForest(as.factor(Species_vec)~., data = dta_Tmp, ntree= 10000)
    prediction_RF_SR[[species]] <- predict(fit_RF_SR[[species]], newdata = Predictor_df_eval, type = "prob")
  }
  
  for(species in 1:100){
    models_fitting$Random_Forest_SR[[trials]][[species]] <- fit_RF_SR[[species]]
  }
  
  #########################################################################################################
  ######################################## Filling 4 D matrix with RF SR ##################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,9,trials] <- prediction_RF_SR[[species]][which(rownames(prediction_RF_SR[[species]]) == plots),2]
      }
  }
  
  #########################################################################################################
  ################################################### GLM5 ################################################
  #########################################################################################################
  
  formula_GLM5 <- vector(length = 100)
  for(species in 1:length(formula_GLM5)){
    formula_GLM5[species] <- paste(sp_names[species], " ~ poly(GS_S_Rad, 2) + poly(GS_p, 2) +  poly(GS_t, 2) + poly(slope, 2) + poly(pH, 2) + poly(sp_richness,2)")
  }
  
  fit_GLM5 <- vector(mode = "list", length = 100)
  prediction_GLM5 <- vector(mode = "list", length = 100)
  for(j in 1:length(fit_GLM5)){
    fit_GLM5[[j]] <- glm(formula_GLM5[j], family = "binomial", data = training_Tmp)
    prediction_GLM5[[j]] <- predict(fit_GLM5[[j]], newdata = evaluating_Tmp, type = "response")
  }
  
  # Save model for Rsquared
  for (species in 1:100) {
    models_fitting$GLM5[[trials]][[species]] <- fit_GLM5[[species]]
  }
  
  #########################################################################################################
  ########################################## Filling 4 D matrix with GLM5 #################################
  #########################################################################################################
  
  for (plots in 1:911) {
    if(plots %in% evaluating_plots)
      for (species in 1:100) {
        forDMatrix[plots,species,10,trials] <- prediction_GLM5[[species]][which(names(prediction_GLM5[[species]]) == plots)]
      }
  }
  
  # The end
}


#########################################################################################################
################################### Creating Three D Matrix for averaging ###############################
#########################################################################################################

dimensions <- c(911, 100, 11)
dinames <- list("Plots", "Species", "Models")
threeDMatrix <- array(data = NA, dim = dimensions)

#Setting dimensions names
dimnames(threeDMatrix)[[1]] <- as.list(Plants_df_scaled$Plot)
dimnames(threeDMatrix)[[2]] <- as.list(sp_names)
dimnames(threeDMatrix)[[3]] <- c("Observed", "GLM1", "GLM2", "HMSC_uncond", "HMSC_cond", "GLM3", "GLM4", "Random_Forest_corr", "Random_Forest_prev", "Random_Forest_SR", "GLM5")

#################################### Filling observed dim in treeDmatrix ################################

for (plots in 1:911) {
  for (species in 1:100) {
    threeDMatrix[plots,species,1] <- Plants_df_scaled[plots,species + 9]
  }
}

#########################################################################################################
######################################## Averaging 4D matrix in 3D matrix ###############################
#########################################################################################################

threeDMatrix[,,2:11] <- apply(forDMatrix, MARGIN = c(1,2,3), function(x) mean(x, na.rm = TRUE))

#########################################################################################################
######################################### Results TABLE + metrics per Species ###########################
#########################################################################################################

Results <- matrix(nrow = 4, ncol = 6)
colnames(Results) <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
rownames(Results) <- c("Mass TSS", "Max Kappa", "AUC", "Rsquared")

list_species <- list(species = vector(length = 100))
list_metrics <- list(TSS = list_species, Kappa = list_species, AUC = list_species, Rsquared = list_species)
list_models <- list(GLM1 = list_metrics, GLM2 = list_metrics, HMSC_uncond = list_metrics,
                    HMSC_cond = list_metrics, GLM3 = list_metrics, GLM4 = list_metrics)

Rsquared_Tjur <- vector(length = 10)

for(models in 1:6){
  
  if(models %in% c(3,4)){
  
    for (species in 1:100) {
    
      list_models[[models]][[1]][[1]][species] <- ecospat.max.tss(threeDMatrix[,species,models+1], threeDMatrix[,species,1])[[2]][1,2]
      list_models[[models]][[2]][[1]][species] <- ecospat.max.kappa(threeDMatrix[,species,models+1], threeDMatrix[,species,1])[[2]][1,2]
      list_models[[models]][[3]][[1]][species] <- auc(roc(threeDMatrix[,species,models+1], as.factor(threeDMatrix[,species,1])))
      
      for (trials in 1:10) {
        Rsquared_Tjur[trials] <- Rsquared(models_fitting[[models]][[trials]][[1]], averageSp = FALSE)[species]
      }
      list_models[[models]][[4]][[1]][species] <- mean(Rsquared_Tjur)
    }
  }
  
  else if(models %in% c(1,2,5,6)){
    
    for (species in 1:100) {
      
      list_models[[models]][[1]][[1]][species] <- ecospat.max.tss(threeDMatrix[,species,models+1], threeDMatrix[,species,1])[[2]][1,2]
      list_models[[models]][[2]][[1]][species] <- ecospat.max.kappa(threeDMatrix[,species,models+1], threeDMatrix[,species,1])[[2]][1,2]
      list_models[[models]][[3]][[1]][species] <- auc(roc(threeDMatrix[,species,models+1], as.factor(threeDMatrix[,species,1])))
      
      for (trials in 1:10) {
        Rsquared_Tjur[trials] <- PseudoR2(models_fitting[[models]][[trials]][[species]], which = "Tjur")
      }
      
      list_models[[models]][[4]][[1]][species] <- mean(Rsquared_Tjur)
      
    }
    
  }
      
  Results[1,models] <- mean(as.numeric(unlist(list_models[[models]]$TSS)))
  Results[2,models] <- mean(as.numeric(unlist(list_models[[models]]$Kappa)))
  Results[3,models] <- mean(as.numeric(unlist(list_models[[models]]$AUC)))
  Results[4,models] <- mean(as.numeric(unlist(list_models[[models]]$Rsquared)))
  
}

#########################################################################################################
########### Statistical Tests to see if models differents between SPECIES PREPARING THE DATA ############
#########################################################################################################

# Creating empty matrices
Stat_analysis_df <- list(TSS = matrix(nrow = 100, ncol = 7), Kappa = matrix(nrow = 100, ncol = 7), AUC = matrix(nrow = 100, ncol = 7), Rsquared = matrix(nrow = 100, ncol = 7))

# Setting columns names and filling the first column with sp_names
for (metrics in 1:4) {
  colnames(Stat_analysis_df[[metrics]]) <- c("Species", "GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
  Stat_analysis_df[[metrics]][,1] <- sp_names 
}

# Filling the matrices with the metrics data
for (metrics in 1:4) {
  for (models in 1:6) {
    Stat_analysis_df[[metrics]][,models+1] <- as.numeric(unlist(list_models[[models]][[metrics]]$species))
  }
}

# Checking metrics and models' metrics distributions
for (metrics in 1:4){
  for (columns in 1:6) {
    Stat_analysis_df[[metrics]][,columns + 1] <- as.numeric(as.character(Stat_analysis_df[[metrics]][, columns + 1]))
    hist(as.numeric(as.character(Stat_analysis_df[[metrics]][,columns + 1])))
    print(shapiro.test(as.numeric(as.character(Stat_analysis_df[[metrics]][,columns + 1]))))
  }
}
# Not normal at all !

# Creating df for lmer (one column with Species, one woth Metrics and one with Models)
lmer_analyses_df <- list(TSS = data.frame(), Kappa = data.frame(), AUC = data.frame(), Rsquared = data.frame())
for (metrics in 1:4){
  for (models in 1:6) {
    if(models == 1)
      lmer_analyses_df[[metrics]] <- cbind(Stat_analysis_df[[metrics]][,1], Stat_analysis_df[[metrics]][,models+1], rep(colnames(Stat_analysis_df[[metrics]])[models + 1], 100))
    else if (models > 1)
      lmer_analyses_df[[metrics]] <- rbind(lmer_analyses_df[[metrics]],cbind(Stat_analysis_df[[metrics]][,1], Stat_analysis_df[[metrics]][,models+1], rep(colnames(Stat_analysis_df[[metrics]])[models + 1], 100)))
  }
  colnames(lmer_analyses_df[[metrics]]) <- c("Species", "Metrics", "Model")
  lmer_analyses_df[[metrics]] <- as.data.frame(lmer_analyses_df[[metrics]])
  lmer_analyses_df[[metrics]]$Model <- factor(lmer_analyses_df[[metrics]]$Model, levels = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"))
  lmer_analyses_df[[metrics]][,2] <- as.numeric(as.character(lmer_analyses_df[[metrics]][,2]))
}

#########################################################################################################
####################################### REPEATED MEASURES ANOVA SPECIES #################################
#########################################################################################################

lm_analyses <- vector(mode = "list", length = 4)
Anova_analyses <- vector(mode = "list", length = 4)

for (metrics in 1:4) {
  lm_analyses[[metrics]] <- lm(Metrics ~ Model + Species, data = lmer_analyses_df[[metrics]])
  Anova_analyses[[metrics]] <- anova(lm_analyses[[metrics]])
}

aov_analyses <- vector(mode = "list", length = 4)
Tukey_analyses <- vector(mode = "list", length = 4)
for (metrics in 1:4) {
  aov_analyses[[metrics]] <- aov(Metrics~Model + Species, data = lmer_analyses_df[[metrics]])
  Tukey_analyses[[metrics]] <- TukeyHSD(aov_analyses[[metrics]])
}

Good_p_values <- vector(mode = "list", length = 4)
for (metrics in 1:4) {
  Good_p_values[[metrics]] <- Tukey_analyses[[metrics]]$Model[,4][Tukey_analyses[[metrics]]$Model[,4] < 0.05]
}

RP_values <- Tukey_analyses[[4]]$Model[,4]
RP_values_corr <- RP_values*4

TSSP_values <- Tukey_analyses[[1]]$Model[,4]
TSSP_values_corr <- TSSP_values*4

KappaP_values <- Tukey_analyses[[2]]$Model[,4]
KappaP_values_corr <- KappaP_values*4

AUCP_values <- Tukey_analyses[[3]]$Model[,4]
AUCP_values_corr <- AUCP_values*4

for (metrics in c(RP_values_corr,TSSP_values_corr,KappaP_values_corr,AUCP_values_corr)) {
  for (values in 1:length(metrics)) {
    if(metrics[values] > 1)
      metrics[values] <- 1
  }
}

#########################################################################################################
############################################## Boxplots metrics for SPECIES #############################
#########################################################################################################

Rp <- ggplot(data = lmer_analyses_df[[4]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "Tjur pseudo-Rsquared", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

TSSp <- ggplot(data = lmer_analyses_df[[1]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "Maximum TSS", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

Kappap <- ggplot(data = lmer_analyses_df[[2]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "Maximum Kappa", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

AUCp <- ggplot(data = lmer_analyses_df[[3]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "AUC", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

grid.arrange(Rp, TSSp, Kappap, AUCp, nrow = 2, ncol = 2)

#########################################################################################################
######################################### Results TABLE + metrics per Plots #############################
#########################################################################################################

Results_plots <- matrix(nrow = 2, ncol = 6)
colnames(Results_plots) <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
rownames(Results_plots) <- c("Community AUC", "Sorensen Index")

list_plots <- list(plots = vector(length = 911))
list_metrics_plots <- list(Commun_AUC = list_plots, Sorensen = list_plots)
list_models_plots <- list(GLM1 = list_metrics_plots, GLM2 = list_metrics_plots, HMSC_uncond = list_metrics_plots,
                          HMSC_cond = list_metrics_plots, GLM3 = list_metrics_plots, GLM4 = list_metrics_plots)

#The MaxSorensen function

MaxSorensen <- function(data){
  obs.data <- as.numeric(data[1:(length(data)/2)])
  pred.data <- as.numeric(data[((length(data)/2)+1):length(data)])
  temp.Sorensen <- rep(NA,101)
  th <- seq(0,1,0.01)
  for(i in 1:101){
    pred.temp <- pred.data
    pred.temp[pred.temp>=th[i]] <- 1
    pred.temp[pred.temp<th[i]] <- 0
    errors <- 2*pred.temp+obs.data
    a <- length(which(errors == 3)) #True Positives
    b <- length(which(errors == 2)) #False Positives
    c <- length(which(errors == 1)) #False Negatives
    if(a==0 & b==0 & c==0){
      Sorensen <- 1
    }else{
      Sorensen <- round((2 * a)/(2 * a + b + c), digits=3)
    }
    temp.Sorensen[i] <- Sorensen
  }
  return(max(temp.Sorensen))
}


for (models in 1:6) {
  for (plots in 1:911) {
    
    if(mean(threeDMatrix[plots,,1]) == 0) {
      list_models_plots[[models]][[1]][[1]][plots] <- NA }
    else {
      list_models_plots[[models]][[1]][[1]][plots] <- auc(roc(threeDMatrix[plots,,models+1], as.factor(threeDMatrix[plots,,1]))) }
    
    list_models_plots[[models]][[2]][[1]][plots] <- MaxSorensen(as.numeric(as.vector(append(threeDMatrix[plots,,1], threeDMatrix[plots,,models+1]))))
    
  }
  
  Results_plots[1,models] <- mean(as.numeric(unlist(list_models_plots[[models]]$Commun_AUC)), na.rm = TRUE)
  Results_plots[2,models] <- mean(as.numeric(unlist(list_models_plots[[models]]$Sorensen)))
  
}

#########################################################################################################
########### Statistical Tests to see if models differents between PLOTS PREPARING THE DATA ##############
#########################################################################################################

# Creating empty matrices
Stat_analysis_df_plots <- list(Commun_AUC = matrix(nrow = 911, ncol = 7), Max_Sorensen = matrix(nrow = 911, ncol = 7))

# Setting columns names and filling the first column with sp_names
for (metrics in 1:2) {
  colnames(Stat_analysis_df_plots[[metrics]]) <- c("Plots", "GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
  Stat_analysis_df_plots[[metrics]][,1] <- Plants_df_scaled$Plot 
}

# Filling the matrices with the metrics data
for (metrics in 1:2) {
  for (models in 1:6) {
    Stat_analysis_df_plots[[metrics]][,models+1] <- as.numeric(unlist(list_models_plots[[models]][[metrics]]$plots))
  }
}

# Checking metrics and models' metrics distributions
for (metrics in 1:2){
  for (columns in 1:6) {
    Stat_analysis_df_plots[[metrics]][,columns + 1] <- as.numeric(as.character(Stat_analysis_df_plots[[metrics]][, columns + 1]))
    hist(as.numeric(as.character(Stat_analysis_df_plots[[metrics]][,columns + 1])))
    print(shapiro.test(as.numeric(as.character(Stat_analysis_df_plots[[metrics]][,columns + 1]))))
  }
}
# Not normal at all !

# Creating df for lmer (one column with Species, one woth Metrics and one with Models)
lmer_analyses_df_plots <- list(CommunAUC = data.frame(), Max_Sorensen = data.frame())
for (metrics in 1:2){
  for (models in 1:6) {
    if(models == 1)
      lmer_analyses_df_plots[[metrics]] <- cbind(Stat_analysis_df_plots[[metrics]][,1], Stat_analysis_df_plots[[metrics]][,models+1], rep(colnames(Stat_analysis_df_plots[[metrics]])[models + 1], 911))
    else if (models > 1)
      lmer_analyses_df_plots[[metrics]] <- rbind(lmer_analyses_df_plots[[metrics]],cbind(Stat_analysis_df_plots[[metrics]][,1], Stat_analysis_df_plots[[metrics]][,models+1], rep(colnames(Stat_analysis_df_plots[[metrics]])[models + 1], 911)))
  }
  colnames(lmer_analyses_df_plots[[metrics]]) <- c("Plots", "Metrics", "Model")
  lmer_analyses_df_plots[[metrics]] <- as.data.frame(lmer_analyses_df_plots[[metrics]])
  lmer_analyses_df_plots[[metrics]]$Model <- factor(lmer_analyses_df_plots[[metrics]]$Model, levels = c("HMSC uncond", "HMSC cond", "GLM1",  "GLM3", "GLM2", "GLM4"))
  lmer_analyses_df_plots[[metrics]][,2] <- as.numeric(as.character(lmer_analyses_df_plots[[metrics]][,2]))
}

#########################################################################################################
################################## REPEATED MEASURES ANOVA PLOTS ########################################
#########################################################################################################

lm_analyses_plots <- vector(mode = "list", length = 2)
Anova_analyses_plots <- vector(mode = "list", length = 2)

for (metrics in 1:2) {
  lm_analyses_plots[[metrics]] <- lm(Metrics ~ Model + Plots, data = lmer_analyses_df_plots[[metrics]])
  Anova_analyses_plots[[metrics]] <- anova(lm_analyses_plots[[metrics]])
}

aov_analyses_plots <- vector(mode = "list", length = 2)
Tukey_analyses_plots <- vector(mode = "list", length = 2)
for (metrics in 1:2) {
  aov_analyses_plots[[metrics]] <- aov(Metrics~Model + Plots, data = lmer_analyses_df_plots[[metrics]])
  Tukey_analyses_plots[[metrics]] <- TukeyHSD(aov_analyses_plots[[metrics]])
}

AUCCP_values <- Tukey_analyses_plots[[1]]$Model[,4]
AUCCP_values_corr <- AUCCP_values*2

MaxSorIndP_values <- Tukey_analyses_plots[[2]]$Model[,4]
MaxSorIndP_values_corr <- MaxSorIndP_values*2

#########################################################################################################
############################################## Boxplots metrics for PLOTS ###############################
#########################################################################################################

AUCCp <- ggplot(data = lmer_analyses_df_plots[[1]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "Community AUC", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

MaxSorIndp <- ggplot(data = lmer_analyses_df_plots[[2]], aes(x = Model, y = Metrics, fill = Model)) + geom_boxplot() +
  labs(y = "Maximum Sorensen Index", x = "Models") +
  scale_x_discrete(breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                   labels = c("HMSC Abiotic", "HMSC Biotic", "GLM Abiotic 1", "GLM Biotic 1", "GLM Abiotic 2", "GLM Biotic 2"),
                   name = NULL) +
  scale_fill_manual(values = c("palegreen1", "palegreen3", "dodgerblue", "dodgerblue4", "orange", "orange3"),
                    breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, angle = 23, hjust = 1),
        axis.text.y = element_text(size = 10))

grid.arrange(AUCCp, MaxSorIndp, nrow = 1, ncol = 2)

#########################################################################################################
####################################### Variance partitionning with HMSC models #########################
#########################################################################################################

################################################### HMSC UNCOND SPECIES #################################

# Calculating variance partitioning
variationPartHMSCuncond <- variPart(model_uncond_HMSC, groupX = colnames(model_uncond_HMSC$data$X))
variationPartHMSCuncond <- as.data.frame(variationPartHMSCuncond)
variationPartHMSCuncond <- cbind(seq(1,100), variationPartHMSCuncond)
colnames(variationPartHMSCuncond)[1] <- "Species_nb"
# Creating adapted dataframe !
var_part_long <- gather(variationPartHMSCuncond, "predictors", "var_proportion", GS_S_Rad:slope_sec)
# Fuse 1th and 2th orders predictors
var_part_long$predictors[var_part_long$predictors == "GS_p_sec"] <- "GS_p"
var_part_long$predictors[var_part_long$predictors == "GS_rad_sec"] <- "GS_S_Rad"
var_part_long$predictors[var_part_long$predictors == "GS_t_sec"] <- "GS_t"
var_part_long$predictors[var_part_long$predictors == "pH_sec"] <- "pH"
var_part_long$predictors[var_part_long$predictors == "slope_sec"] <- "slope"

# Plotting
HMSCuncondVarPart <- ggplot(var_part_long, aes(x = as.factor(Species_nb), y = var_proportion, fill = predictors)) +
  geom_col(position = "stack") + labs(x = "Species", y = "Proportion of Variance explained within the model", title = "A  HMSC Abiotic") +
  scale_fill_manual(values = c("#00BFC4", "#E58700", "#FF99CC", "#B983FF", "#CCCC33"), 
                    name = "Predictors", breaks = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope"),
                    labels = c("Precipitations", "Solar Radiations", "Temperature", "pH", "Slope")) +
  scale_x_discrete(name = "Species", breaks = NULL) + guides(fill = FALSE) +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 18, face = "bold"), axis.text = element_text(size = 10), axis.title = element_text(size = 12))

################################################### HMSC COND SPECIES #########################################

variationPartHMSCcond <- variPart(model_cond_HMSC, groupX = colnames(model_cond_HMSC$data$X))
variationPartHMSCcond <- as.data.frame(variationPartHMSCcond)
variationPartHMSCcond <- cbind(seq(1,100), variationPartHMSCcond)
colnames(variationPartHMSCcond)[1] <- "Species_nb"
colnames(variationPartHMSCcond)[13] <- "Random"
colnames(variationPartHMSCcond)[14] <- "Random2"
# Creating adapted dataframe !
var_part_long_cond <- gather(variationPartHMSCcond, "predictors", "var_proportion", GS_S_Rad:Random2)
# Fuse 1th and 2th orders predictors
var_part_long_cond$predictors[var_part_long_cond$predictors == "GS_p_sec"] <- "GS_p"
var_part_long_cond$predictors[var_part_long_cond$predictors == "GS_rad_sec"] <- "GS_S_Rad"
var_part_long_cond$predictors[var_part_long_cond$predictors == "GS_t_sec"] <- "GS_t"
var_part_long_cond$predictors[var_part_long_cond$predictors == "pH_sec"] <- "pH"
var_part_long_cond$predictors[var_part_long_cond$predictors == "slope_sec"] <- "slope"
var_part_long_cond$predictors[var_part_long_cond$predictors == "Random2"] <- "Random"
var_part_long_cond$predictors <- factor(var_part_long_cond$predictors, levels = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"))

# Plotting
HMSCcondVarPart <- ggplot(var_part_long_cond, aes(x = as.factor(Species_nb), y = var_proportion, fill = as.factor(predictors))) +
  geom_col(position = "stack") + labs(x = "Species", y = "Proportion of Variance explained within the model", title = "HMSC Biotic") +
  scale_fill_manual(values = c("#00BFC4", "#E58700", "#FF99CC", "#B983FF", "#00CC66", "#CCCC33"),
                    name = "Predictors", breaks = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"),
                    labels = c("Precipitations", "Solar Radiations", "Temperature", "pH", "Slope", "Latent Variables")) +
  scale_x_discrete(name = "Species", breaks = seq(0,100,5)) + theme_bw() +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
              legend.title=element_text(size=18), 
              legend.text=element_text(size=15))

#########################################################################################################
####################################### SPECIES are better which modelling method #######################
#########################################################################################################

#Create empty vectors for Kappa for each modelling method and for the best model for each species
Kappa <- vector(length = 6)
names(Kappa) <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
Best_method <- vector(length = 100)
Best_method_1 <- vector(length = 100)
names(Best_method) <- sp_names
for(species in 1:100){
  for(models in 1:6){
    Kappa[models] <- list_models[[models]]$Kappa$species[species]
  }
  Best_method_1[species] <- list(names(Kappa[Kappa == max(Kappa)]))
}

for (species in 1:100) {
  if(length(Best_method_1[[species]]) == 1)
    Best_method[species] <- Best_method_1[[species]]
  else
  Best_method[species] <- "GLM4"
}

table_best_method <- table(Best_method)

#########################################################################################################
################################################ Percentages SPECIES ####################################
#########################################################################################################

(37 + 45 + 7)/100 # 89% Species better with BIOTIC MODELS

(37 + 45) / 89 # Within those species better with biotic models, 92% better with GLMs

#########################################################################################################
################################ Similarites between Species better with abiotic ########################
#########################################################################################################

Abiotic_Better <- Best_method[Best_method == "GLM1" | Best_method == "GLM2" | Best_method == "HMSC uncond"]
Biotic_Better <- Best_method[Best_method == "GLM3" | Best_method == "GLM4" | Best_method == "HMSC cond"]

#########################################################################################################
############################################### PREVALENCE ##############################################
#########################################################################################################

Abiotic_Better_Prevalence <- vector(length = 11)
names(Abiotic_Better_Prevalence) <- names(Abiotic_Better)

Biotic_Better_Prevalence <- vector(length = 89)
names(Biotic_Better_Prevalence) <- names(Biotic_Better)

for (names in names(Abiotic_Better)) {
  Plant_Tmp <- Plants_df_scaled[colnames(Plants_df_scaled) == names]
  Abiotic_Better_Prevalence[names(Abiotic_Better_Prevalence) == names] <- length(Plant_Tmp[Plant_Tmp == "1"])
}

for (names in names(Biotic_Better)) {
  Plant_Tmp <- Plants_df_scaled[colnames(Plants_df_scaled) == names]
  Biotic_Better_Prevalence[names(Biotic_Better_Prevalence) == names] <- length(Plant_Tmp[Plant_Tmp == "1"])
}

#########################################################################################################
################################################ LEAF AREA ##############################################
#########################################################################################################
Traits <- read.csv("TraitsVF.csv", sep = ";")

sla <- Traits$AvgOfSLA
names(sla) <- gsub(" " , "_", Traits$spnames, fixed = TRUE)

sla <- sla[names(sla) %in% sp_names]
missing <- sp_names[! sp_names %in% names(sla)]

Abiotic_better_LA <- vector(length = 11)
names(Abiotic_better_LA) <- names(Abiotic_Better)

for (names in names(Abiotic_Better)) {
  if(names %in% names(sla))
    Abiotic_better_LA[names(Abiotic_better_LA) == names] <- sla[names(sla) == names]
}

Biotic_Better_LA <- vector(length = 89)
names(Biotic_Better_LA) <- names(Biotic_Better)

for (names in names(Biotic_Better)) {
  if(names %in% names(sla))
    Biotic_Better_LA[names(Biotic_Better_LA) == names] <- sla[names(sla) == names]
}

#########################################################################################################
###################################### BOXPLOTS SPECIES better modelling method #########################
#########################################################################################################
Modeling_Better_species <- c(rep("Biotic",89), rep("Abiotic",11))
Modeling_Better_species <- cbind(Modeling_Better_species, append(Biotic_Better, Abiotic_Better),
                                 append(Biotic_Better_Prevalence, Abiotic_Better_Prevalence),
                                 append(Biotic_Better_LA, Abiotic_better_LA))
Modeling_Better_species <- as.data.frame(Modeling_Better_species)
colnames(Modeling_Better_species) <- c("Best_Method", "Best_Model", "Prevalence", "Leaf_Area")

Prev <- ggplot(data = Modeling_Better_species, aes(x = Best_Method)) +
  geom_boxplot(aes(y = as.numeric(as.character(Prevalence))), fill = "lightgrey") + labs(x = NULL, y = "Prevalence") +
  scale_x_discrete(breaks = c("Biotic", "Abiotic"),
                   labels = c("Species better with Biotic models", "Species better with Abiotic models")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14))

wilcox.test(as.numeric(as.character(Modeling_Better_species$Prevalence[Modeling_Better_species$Best_Method == "Biotic"])),
            as.numeric(as.character(Modeling_Better_species$Prevalence[Modeling_Better_species$Best_Method == "Abiotic"])))

LA <- ggplot(data = Modeling_Better_species, aes(x = Best_Method)) +
  geom_boxplot(aes(y = as.numeric(as.character(Leaf_Area))), fill = "lightgrey") + labs(x = NULL, y = "Specific Leaf Area") +
  scale_x_discrete(breaks = c("Biotic", "Abiotic"),
                   labels = c("Species better with Biotic models", "Species better with Abiotic models")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14))

wilcox.test(as.numeric(as.character(Modeling_Better_species$Leaf_Area[Modeling_Better_species$Best_Method == "Biotic"])),
            as.numeric(as.character(Modeling_Better_species$Leaf_Area[Modeling_Better_species$Best_Method == "Abiotic"])))

grid.arrange(Prev, LA, nrow = 1, ncol = 2)

#########################################################################################################
####################################### PLOTS are better which modelling method #########################
#########################################################################################################

#Create empty vectors for Kappa for each modelling method and for the best model for each species
Sorensen <- vector(length = 6)
names(Sorensen) <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
Best_method_plots_1 <- vector(length = 911)
names(Best_method_plots_1) <- seq(1,911)
Best_method_plots <- vector(length = 911)
names(Best_method_plots) <- seq(1,911)
for(plots in 1:911){
  for(models in 1:6){
    Sorensen[models] <- list_models_plots[[models]]$Sorensen$plots[plots]
  }
  Best_method_plots_1[plots] <- list(names(Sorensen[Sorensen == max(Sorensen)]))
}

for (plots in 1:911) {
  if(length(Best_method_plots_1[[plots]]) == 1){
    Best_method_plots[plots] <- Best_method_plots_1[[plots]]}
  else if(length(Best_method_plots_1[[plots]]) > 1){
    if((TRUE %in% (Best_method_plots_1[[plots]] %in% c("GLM1", "GLM2", "HMSC uncond"))) & (TRUE %in% (Best_method_plots_1[[plots]] %in% c("GLM3", "GLM4", "HMSC cond")))){
      Best_method_plots[plots] <- "FALSE"}
    else if((TRUE %in% (Best_method_plots_1[[plots]] %in% c("GLM1", "GLM2", "HMSC uncond"))) & (FALSE %in% (Best_method_plots_1[[plots]] %in% c("GLM3", "GLM4", "HMSC cond")))){
      Best_method_plots[plots] <- "GLM1"}
    else if((FALSE %in% (Best_method_plots_1[[plots]] %in% c("GLM1", "GLM2", "HMSC uncond"))) & (TRUE %in% (Best_method_plots_1[[plots]] %in% c("GLM3", "GLM4", "HMSC cond")))){
      Best_method_plots[plots] <- "GLM4"}
  }
}

table_best_method_plots <- table(Best_method_plots)

#########################################################################################################
################################################ Percentages PLOTS ####################################
#########################################################################################################

222/911 # 24% Plots with several modelling methods

(213 + 187 + 74)/(911-222) # 69% Plots better with BIOTIC MODELS within the taken into account plots

(213 + 187) / 474 # Within those Plots better with biotic models, 84% better with GLMs

#########################################################################################################
################################ Similarites between PLOTS better with abiotic ##########################
#########################################################################################################

Abiotic_Better_plots <- Best_method_plots[Best_method_plots == "GLM1" | Best_method_plots == "GLM2" | Best_method_plots == "HMSC uncond"]
Biotic_Better_plots <- Best_method_plots[Best_method_plots == "GLM3" | Best_method_plots == "GLM4" | Best_method_plots == "HMSC cond"]


#########################################################################################################
################################################### ELEVATION ###########################################
#########################################################################################################

Abiotic_Better_plots_elevation <- vector(length = 215)
names(Abiotic_Better_plots_elevation) <- names(Abiotic_Better_plots)

for (names in names(Abiotic_Better_plots)) {
  Abiotic_Better_plots_elevation[names(Abiotic_Better_plots_elevation) == names] <- Plants_df_scaled$DEM_vec[Plants_df_scaled$Plot == names]
}

Biotic_Better_plots_elevation <- vector(length = 474)
names(Biotic_Better_plots_elevation) <- names(Biotic_Better_plots)

for (names in names(Biotic_Better_plots)) {
  Biotic_Better_plots_elevation[names(Biotic_Better_plots_elevation) == names] <- Plants_df_scaled$DEM_vec[Plants_df_scaled$Plot == names]
}

#########################################################################################################
############################################ SPECIES RICHNESS ###########################################
#########################################################################################################

Abiotic_Better_plots_SR <- vector(length = 215)
names(Abiotic_Better_plots_SR) <- names(Abiotic_Better_plots)

for (names in names(Abiotic_Better_plots)) {
  Abiotic_Better_plots_SR[names(Abiotic_Better_plots_SR) == names] <- sum(Plants_df_scaled[Plants_df_scaled$Plot == names,10:109])
}

Biotic_Better_plots_SR <- vector(length = 474)
names(Biotic_Better_plots_SR) <- names(Biotic_Better_plots)

for (names in names(Biotic_Better_plots)) {
  Biotic_Better_plots_SR[names(Biotic_Better_plots_SR) == names] <- sum(Plants_df_scaled[Plants_df_scaled$Plot == names, 10:109])
}

#########################################################################################################
################################## BOXPLOTS PLOTS better modelling methods ##############################
#########################################################################################################

Modeling_Better_plots <- c(rep("Biotic",474), rep("Abiotic",215))
Modeling_Better_plots <- cbind(Modeling_Better_plots, append(Biotic_Better_plots, Abiotic_Better_plots),
                                 append(Biotic_Better_plots_elevation, Abiotic_Better_plots_elevation),
                                 append(Biotic_Better_plots_SR, Abiotic_Better_plots_SR))

Modeling_Better_plots <- as.data.frame(Modeling_Better_plots)
colnames(Modeling_Better_plots) <- c("Best_Method", "Best_Model", "Elevation", "Species_Richness")

Elev_plots <- ggplot(data = Modeling_Better_plots, aes(x = Best_Method)) +
  geom_boxplot(aes(y = as.numeric(as.character(Elevation))), fill = "lightgrey") + labs(x = NULL, y = "Elevation") +
  scale_x_discrete(breaks = c("Biotic", "Abiotic"),
                   labels = c("Species better with Biotic models", "Species better with Abiotic models")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14))

wilcox.test(as.numeric(as.character(Modeling_Better_plots$Elevation[Modeling_Better_plots$Best_Method == "Biotic"])),
            as.numeric(as.character(Modeling_Better_plots$Elevation[Modeling_Better_plots$Best_Method == "Abiotic"])))

SR <- ggplot(data = Modeling_Better_plots, aes(x = Best_Method)) +
  geom_boxplot(aes(y = as.numeric(as.character(Species_Richness))), fill = "lightgrey") + labs(x = NULL, y = "Species Richness") +
  scale_x_discrete(breaks = c("Biotic", "Abiotic"),
                   labels = c("Species better with Biotic models", "Species better with Abiotic models")) +
  theme_bw() + guides(fill = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14))

wilcox.test(as.numeric(as.character(Modeling_Better_plots$Species_Richness[Modeling_Better_plots$Best_Method == "Biotic"])),
            as.numeric(as.character(Modeling_Better_plots$Species_Richness[Modeling_Better_plots$Best_Method == "Abiotic"])))

grid.arrange(Elev_plots, SR, nrow = 1, ncol = 2)

#########################################################################################################
##################################################### ELEVATION GRAPH  ##################################
#########################################################################################################

# Creating the data frames
Elevation_graph_df <- cbind(unlist(list_models_plots$GLM1$Sorensen), unlist(list_models_plots$GLM2$Sorensen),
                            unlist(list_models_plots$HMSC_uncond$Sorensen), unlist(list_models_plots$HMSC_cond$Sorensen),
                            unlist(list_models_plots$GLM3$Sorensen), unlist(list_models_plots$GLM4$Sorensen))
Elevation_graph_df <- as.data.frame(Elevation_graph_df)
colnames(Elevation_graph_df) <- c("GLM1", "GLM2", "HMSC uncond", "HMSC cond", "GLM3", "GLM4")
Elevation_graph_df <- cbind(Plants_df_scaled$DEM_vec, Plants_df_scaled$Plot, Elevation_graph_df)
colnames(Elevation_graph_df)[1] <- "Elevation"
colnames(Elevation_graph_df)[2] <- "Plots"

Elevation_graph_ordered <- Elevation_graph_df[order(Elevation_graph_df$Elevation),]
Elevation_graph_gathered <- gather(Elevation_graph_ordered, "Models", "Max Sorensen Index", GLM1:GLM4)

min_df <- vector(length = 911)
for (plots in 1:911) {
  min_df[plots] <- min(Elevation_graph_ordered[plots,3:8])
}

max_df <- vector(length = 911)
for (plots in 1:911) {
  max_df[plots] <- max(Elevation_graph_ordered[plots,3:8])
}

Elevation_graph_gathered_mm <- cbind(Elevation_graph_gathered, append(min_df, min_df), append(max_df, max_df))
colnames(Elevation_graph_gathered_mm)[5] <- "Min"
colnames(Elevation_graph_gathered_mm)[6] <- "Max"
Elevation_graph_gathered_mm$Models <- factor(Elevation_graph_gathered_mm$Models, levels = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"))
Elevation_graph_gathered$Models <- factor(Elevation_graph_gathered$Models, levels = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"))
dta <- Plants_df_scaled[,c(9,110)]
colnames(dta)[1] <- "Elevation"
Elevation_graph_gathered_mm <- merge(Elevation_graph_gathered_mm, dta, by = "Elevation")

elev_plot <- ggplot(data = Elevation_graph_gathered_mm, aes(x = Elevation, y = as.numeric(as.character(`Max Sorensen Index`)), colour = Models)) +
  labs(x = "Elevation [m]", y = "Maximum Sorensen Index") + geom_point(size = 5) +
  geom_segment(aes(x = Elevation, xend = Elevation, y = Min, yend = Max), col = "black") +
  scale_color_manual(values = c("GLM1" = "dodgerblue", "GLM2" = "orange", "HMSC uncond" = "palegreen1", "HMSC cond" = "palegreen3", "GLM3" = "dodgerblue4", "GLM4" = "orange3"),
                     breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                     labels = c("HMSC Abiotic", "HMSC Biotic", "Abiotic GLM1", "Biotic GLM1", "Abiotic GLM2", "Biotic GLM2")) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text=element_text(size=20), axis.title=element_text(size=24))

#########################################################################################################
######################################## ELEVATION GRAPH TESTS ##########################################
#########################################################################################################

anova(lm(`Max Sorensen Index`~ Models + Elevation, data = Elevation_graph_gathered))

#########################################################################################################
################################## LOCAL POLYNOMIAL REGRESSION + SR #####################################
#########################################################################################################
Plant_df_scaled_ordered <- Plants_df_scaled[order(Plants_df_scaled$DEM_vec),]
Plant_df_scaled_ordered$Plot <- factor(Plant_df_scaled_ordered$Plot, levels = Plant_df_scaled_ordered$Plot)

pmod <- ggplot(data = Elevation_graph_gathered, aes(x = Elevation, y = `Max Sorensen Index`, colour = Models)) +
  geom_smooth(method = "loess", se = TRUE) + labs(x = "Elevation [m]", y = "Maximum Sorensen Index") + ylim(c(0,1)) +
  scale_color_manual(values = c("GLM1" = "dodgerblue", "GLM2" = "orange", "HMSC uncond" = "palegreen1", "HMSC cond" = "palegreen3", "GLM3" = "dodgerblue4", "GLM4" = "orange3"),
                     breaks = c("HMSC uncond", "HMSC cond", "GLM1", "GLM3", "GLM2", "GLM4"),
                     labels = c("HMSC Abiotic", "HMSC Biotic", "Abiotic GLM1", "Biotic GLM1", "Abiotic GLM2", "Biotic GLM2")) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.text=element_text(size=16),
        axis.title=element_text(size=18))

pSR <- ggplot(data = Elevation_graph_gathered_mm, aes(x = Elevation, y = sp_richness)) +
  geom_smooth(data = Plant_df_scaled_ordered, aes(x = DEM_vec, y = sp_richness),
              method = "loess", col = "black", span = 0.6, se = TRUE) + 
  labs(x = "Elevation [m]", y = "Species Richness") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.text=element_text(size=16),
        axis.title=element_text(size=18))

g1 <- ggplot_gtable(ggplot_build(pmod))
g2 <- ggplot_gtable(ggplot_build(pSR))

pp <- c(subset(g1$layout, name == "panel", se = t:r))

g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]

ax$widths <- rev(ax$widths)

g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

ggsave("plot.pdf", g, width = 10, height = 6)

#########################################################################################################
################################# VARIANCE PARTIONNING SEPARATED ELEVATION ##############################
#########################################################################################################

#################################### Creating low, mid and high elev ####################################
Plants_df_scaled_Low <- Plant_df_scaled_ordered[1:304,]
Plants_df_scaled_Mid <- Plant_df_scaled_ordered[305:607,]
Plants_df_scaled_High <- Plant_df_scaled_ordered[608:911,]

#########################################################################################################
############################################### LOW ELEV ################################################

# Building appropriate matrices
Low_elev_plant <- Plants_df_scaled_Low[,10:109]
Low_elev_plant <- as.matrix(Low_elev_plant)
Low_elev_var <- Plants_df_scaled_Low[,4:8]
second_term_mat_eval <- matrix(ncol = 5, nrow = nrow(Plants_df_scaled_Low))
second_term_mat_eval <- as.data.frame(second_term_mat_eval)
second_term_mat_eval <- sapply(Low_elev_var, function(x) x^2)
colnames(second_term_mat_eval) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
Low_elev_var <- cbind(Low_elev_var, second_term_mat_eval)
Low_elev_var <- as.matrix(Low_elev_var)

# Ramdom matrix
Random_factor_Low <- rep(1,nrow(Plants_df_scaled_Low))
Random_factor_Low <- as.factor(Random_factor_Low)

# Formatting data
Form_low_elev <- as.HMSCdata(Y = Low_elev_plant, X = Low_elev_var, Random = Random_factor_Low, interceptX = TRUE, scaleX = FALSE)

# Fitting the model
model_cond_Low <- hmsc(Form_low_elev, family = "probit", niter = 10000, nburn = 5000, thin = 50)

# VARIATION GRAPH
variationPartHMSCcond_Low <- variPart(model_cond_Low, groupX = colnames(model_cond_Low$data$X))
variationPartHMSCcond_Low <- as.data.frame(variationPartHMSCcond_Low)
variationPartHMSCcond_Low <- cbind(seq(1,100), variationPartHMSCcond_Low)
colnames(variationPartHMSCcond_Low)[1] <- "Species_nb"
colnames(variationPartHMSCcond_Low)[13] <- "Random"
# Creating adapted dataframe !
var_part_long_cond_Low <- gather(variationPartHMSCcond_Low, "predictors", "var_proportion", GS_S_Rad:Random)
# Fuse 1th and 2th orders predictors
var_part_long_cond_Low$predictors[var_part_long_cond_Low$predictors == "GS_p_sec"] <- "GS_p"
var_part_long_cond_Low$predictors[var_part_long_cond_Low$predictors == "GS_rad_sec"] <- "GS_S_Rad"
var_part_long_cond_Low$predictors[var_part_long_cond_Low$predictors == "GS_t_sec"] <- "GS_t"
var_part_long_cond_Low$predictors[var_part_long_cond_Low$predictors == "pH_sec"] <- "pH"
var_part_long_cond_Low$predictors[var_part_long_cond_Low$predictors == "slope_sec"] <- "slope"
var_part_long_cond_Low$predictors <- factor(var_part_long_cond_Low$predictors, levels = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"))

# Plotting
ggplot(var_part_long_cond_Low, aes(x = as.factor(Species_nb), y = var_proportion, fill = as.factor(predictors))) +
  geom_col(position = "stack") + labs(x = "Species", y = "Proportion of Variance explained within the model", title = "Low Elevation (418.8m ; 1479.4m)") +
  scale_fill_manual(values = c("#00BFC4", "#E58700", "#FF99CC", "#B983FF", "#CCCC33", "#00CC66"),
                    name = "Predictors", breaks = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"),
                    labels = c("Precipitations", "Solar Radiations", "Temperature", "pH", "Slope", "Random")) +
  scale_x_discrete(name = "Species", breaks = seq(0,100,10)) +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 12, face = "bold"), legend.position = "none",
        legend.title=element_text(size=18), legend.text=element_text(size=15))

#########################################################################################################
############################################### MID ELEV ################################################

# Building appropriate matrices
Mid_elev_plant <- Plants_df_scaled_Mid[,10:109]
Mid_elev_plant <- as.matrix(Mid_elev_plant)
Mid_elev_var <- Plants_df_scaled_Mid[,4:8]
second_term_mat_eval <- matrix(ncol = 5, nrow = nrow(Plants_df_scaled_Mid))
second_term_mat_eval <- as.data.frame(second_term_mat_eval)
second_term_mat_eval <- sapply(Mid_elev_var, function(x) x^2)
colnames(second_term_mat_eval) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
Mid_elev_var <- cbind(Mid_elev_var, second_term_mat_eval)
Mid_elev_var <- as.matrix(Mid_elev_var)

# Ramdom matrix
Random_factor_Mid <- rep(1,nrow(Plants_df_scaled_Mid))
Random_factor_Mid <- as.factor(Random_factor_Mid)

# Formatting data
Form_Mid_elev <- as.HMSCdata(Y = Mid_elev_plant, X = Mid_elev_var, Random = Random_factor_Mid, interceptX = TRUE, scaleX = FALSE)

# Fitting the model
model_cond_Mid <- hmsc(Form_Mid_elev, family = "probit", niter = 10000, nburn = 5000, thin = 50)

# VARIATION GRAPH
variationPartHMSCcond_Mid <- variPart(model_cond_Mid, groupX = colnames(model_cond_Mid$data$X))
variationPartHMSCcond_Mid <- as.data.frame(variationPartHMSCcond_Mid)
variationPartHMSCcond_Mid <- cbind(seq(1,100), variationPartHMSCcond_Mid)
colnames(variationPartHMSCcond_Mid)[1] <- "Species_nb"
colnames(variationPartHMSCcond_Mid)[13] <- "Random"
# Creating adapted dataframe !
var_part_long_cond_Mid <- gather(variationPartHMSCcond_Mid, "predictors", "var_proportion", GS_S_Rad:Random)
# Fuse 1th and 2th orders predictors
var_part_long_cond_Mid$predictors[var_part_long_cond_Mid$predictors == "GS_p_sec"] <- "GS_p"
var_part_long_cond_Mid$predictors[var_part_long_cond_Mid$predictors == "GS_rad_sec"] <- "GS_S_Rad"
var_part_long_cond_Mid$predictors[var_part_long_cond_Mid$predictors == "GS_t_sec"] <- "GS_t"
var_part_long_cond_Mid$predictors[var_part_long_cond_Mid$predictors == "pH_sec"] <- "pH"
var_part_long_cond_Mid$predictors[var_part_long_cond_Mid$predictors == "slope_sec"] <- "slope"
var_part_long_cond_Mid$predictors <- factor(var_part_long_cond_Mid$predictors, levels = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"))

# Plotting
ggplot(var_part_long_cond_Mid, aes(x = as.factor(Species_nb), y = var_proportion, fill = as.factor(predictors))) +
  geom_col(position = "stack") + labs(x = "Species", y = "Proportion of Variance explained within the model", title = "Mid Elevation (1479.4m ; 2040.5m)") +
  scale_fill_manual(values = c("#00BFC4", "#E58700", "#FF99CC", "#B983FF", "#CCCC33", "#00CC66"),
                    name = "Predictors", breaks = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"),
                    labels = c("Precipitations", "Solar Radiations", "Temperature", "pH", "Slope", "Random")) +
  scale_x_discrete(name = "Species", breaks = seq(0,100,10)) +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 12, face = "bold"), legend.position = "none",
        legend.title=element_text(size=18), legend.text=element_text(size=15))

#########################################################################################################
############################################### HIGH ELEV ################################################

# Building appropriate matrices
High_elev_plant <- Plants_df_scaled_High[,10:109]
High_elev_plant <- as.matrix(High_elev_plant)
High_elev_var <- Plants_df_scaled_High[,4:8]
second_term_mat_eval <- matrix(ncol = 5, nrow = nrow(Plants_df_scaled_High))
second_term_mat_eval <- as.data.frame(second_term_mat_eval)
second_term_mat_eval <- sapply(High_elev_var, function(x) x^2)
colnames(second_term_mat_eval) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
High_elev_var <- cbind(High_elev_var, second_term_mat_eval)
High_elev_var <- as.matrix(High_elev_var)

# Ramdom matrix
Random_factor_High <- rep(1,nrow(Plants_df_scaled_High))
Random_factor_High <- as.factor(Random_factor_High)

# Formatting data
Form_High_elev <- as.HMSCdata(Y = High_elev_plant, X = High_elev_var, Random = Random_factor_High, interceptX = TRUE, scaleX = FALSE)

# Fitting the model
model_cond_High <- hmsc(Form_High_elev, family = "probit", niter = 10000, nburn = 5000, thin = 50)

# VARIATION GRAPH
variationPartHMSCcond_High <- variPart(model_cond_High, groupX = colnames(model_cond_High$data$X))
variationPartHMSCcond_High <- as.data.frame(variationPartHMSCcond_High)
variationPartHMSCcond_High <- cbind(seq(1,100), variationPartHMSCcond_High)
colnames(variationPartHMSCcond_High)[1] <- "Species_nb"
colnames(variationPartHMSCcond_High)[13] <- "Random"
# Creating adapted dataframe !
var_part_long_cond_High <- gather(variationPartHMSCcond_High, "predictors", "var_proportion", GS_S_Rad:Random)
# Fuse 1th and 2th orders predictors
var_part_long_cond_High$predictors[var_part_long_cond_High$predictors == "GS_p_sec"] <- "GS_p"
var_part_long_cond_High$predictors[var_part_long_cond_High$predictors == "GS_rad_sec"] <- "GS_S_Rad"
var_part_long_cond_High$predictors[var_part_long_cond_High$predictors == "GS_t_sec"] <- "GS_t"
var_part_long_cond_High$predictors[var_part_long_cond_High$predictors == "pH_sec"] <- "pH"
var_part_long_cond_High$predictors[var_part_long_cond_High$predictors == "slope_sec"] <- "slope"
var_part_long_cond_High$predictors <- factor(var_part_long_cond_High$predictors, levels = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"))

# Plotting
ggplot(var_part_long_cond_High, aes(x = as.factor(Species_nb), y = var_proportion, fill = as.factor(predictors))) +
  geom_col(position = "stack") + labs(x = "Species", y = "Proportion of Variance explained within the model", title = "High Elevation (2040.7m ; 3101.4m)") +
  scale_fill_manual(values = c("#00BFC4", "#E58700", "#FF99CC", "#B983FF", "#CCCC33", "#00CC66"),
                    name = "Predictors", breaks = c("GS_p", "GS_S_Rad", "GS_t", "pH", "slope", "Random"),
                    labels = c("Precipitations", "Solar Radiations", "Temperature", "pH", "Slope", "Biotic Variables (Latent Variables)")) +
  scale_x_discrete(name = "Species", breaks = seq(0,100,5)) +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=18), legend.text=element_text(size=15))

#########################################################################################################
############################# Prop test VARIANCE PARTITIONING ALONG ELEVATION ###########################
#########################################################################################################
Ratio_Low <- vector(length = 100)
for(species in 1:100) {
  abiotic <- sum(var_part_long_cond_Low[species,3:12])
  random <- var_part_long_cond_Low[species,13]
  Ratio_Low[species] <- random/abiotic
}

Ratio_Mid <- vector(length = 100)
for(species in 1:100) {
  abiotic <- sum(var_part_long_cond_Mid[species,3:12])
  random <- var_part_long_cond_Mid[species,13]
  Ratio_Mid[species] <- random/abiotic
}

Ratio_High <- vector(length = 100)
for(species in 1:100) {
  abiotic <- sum(var_part_long_cond_High[species,3:12])
  random <- var_part_long_cond_High[species,13]
  Ratio_High[species] <- random/abiotic
}

Ratio_df <- cbind(seq(1,100),Ratio_Low,Ratio_Mid,Ratio_High)
colnames(Ratio_df) <- c("Species", "Low", "Mid", "High")
Ratio_df <- as.data.frame(Ratio_df)

Ratio_df_gathered <- gather(Ratio_df, key = "Elevation", value = "Ratio", Low:High)
Ratio_df_gathered$Species <- as.factor(Ratio_df_gathered$Species)

model_elev <- lm(Ratio~Elevation + Species, data = Ratio_df_gathered)
anova(model_elev)
model_elev_2 <- aov(Ratio~Elevation + Species, data = Ratio_df_gathered)
Tukey_elev <- TukeyHSD(model_elev_2)

#########################################################################################################
############################################## ADDITIONAL ANALYSES ######################################
#########################################################################################################

#########################################################################################################
#####################################  Random Forest 2 ##################################################
#########################################################################################################

Results_RF <- vector(length = 4)
names(Results_RF) <- c("Mass.TSS", "Max.Kappa", "AUC", "Rsquared")

list_species <- list(species = vector(length = 100))
RF_Comparison <- matrix(ncol = 6, nrow = 100)
colnames(RF_Comparison) <- c("Species", "TSS", "Kappa", "AUC", "Rsquared", "Model")
RF_Comparison[,1] <- seq(1,100)
RF_Comparison[,6] <- rep("Random.Forest.prev",100)

Rsquared_Tjur_RF <- vector(length = 10)

for (species in 1:100) {
  
  RF_Comparison[species,2] <- ecospat.max.tss(threeDMatrix[,species,9], threeDMatrix[,species,1])[[2]][1,2]
  RF_Comparison[species,3] <- ecospat.max.kappa(threeDMatrix[,species,9], threeDMatrix[,species,1])[[2]][1,2]
  RF_Comparison[species,4] <- auc(roc(threeDMatrix[,species,9], as.factor(threeDMatrix[,species,1])))
  
  for (trials in 1:10) {
    Rsquared_Tjur_RF[trials] <- PseudoR2(models_fitting[[8]][[trials]][[species]], which = "Tjur")
  }
  
  RF_Comparison[species,5] <- mean(Rsquared_Tjur_RF, na.rm = TRUE)
}

#########################################################################################################
########################################### Random Forest 1 #############################################
#########################################################################################################

Results_RF_plots <- vector(length = 2)
names(Results_RF_plots) <- c("Commu.AUC", "Max.Sorensen")

list_plots <- list(plots = vector(length = 911))
RF_Comparison_plots <- matrix(ncol = 4, nrow = 911)
colnames(RF_Comparison_plots) <- c("Plots", "CommunAUC", "MaxSorensen", "Model")
RF_Comparison_plots[,1] <- seq(1,911)
RF_Comparison_plots[,4] <- rep("Random.Forest.cor",911)

for (plots in 1:911) {
  
  if(mean(threeDMatrix[plots,,1]) == 0) {
    RF_Comparison_plots[plots,2] <- NA }
  else {
    RF_Comparison_plots[plots,2] <- auc(roc(threeDMatrix[plots,,8], as.factor(threeDMatrix[plots,,1]))) }
  
  RF_Comparison_plots[plots,3] <- MaxSorensen(as.numeric(as.vector(append(threeDMatrix[plots,,1], threeDMatrix[plots,,8]))))
  
}

#########################################################################################################
####################################### RF SPECIES RICHNESS #############################################
#########################################################################################################

Results_RF_SR <- vector(length = 4)
names(Results_RF_SR) <- c("Mass.TSS", "Max.Kappa", "AUC", "Rsquared")

list_species <- list(species = vector(length = 100))
RF_Comparison_SR <- matrix(ncol = 6, nrow = 100)
colnames(RF_Comparison_SR) <- c("Species", "TSS", "Kappa", "AUC", "Rsquared", "Model")
RF_Comparison_SR[,1] <- seq(1,100)
RF_Comparison_SR[,6] <- rep("Random.Forest.SR",100)

Rsquared_Tjur_RF_SR <- vector(length = 10)

for (species in 1:100) {
  
  RF_Comparison_SR[species,2] <- ecospat.max.tss(threeDMatrix[,species,10], threeDMatrix[,species,1])[[2]][1,2]
  RF_Comparison_SR[species,3] <- ecospat.max.kappa(threeDMatrix[,species,10], threeDMatrix[,species,1])[[2]][1,2]
  RF_Comparison_SR[species,4] <- auc(roc(threeDMatrix[,species,10], as.factor(threeDMatrix[,species,1])))
  
  for (trials in 1:10) {
    Rsquared_Tjur_RF[trials] <- PseudoR2(models_fitting[[9]][[trials]][[species]], which = "Tjur")
  }
  
  RF_Comparison_SR[species,5] <- mean(Rsquared_Tjur_RF_SR)
}

Results_RF_plots_SR <- vector(length = 2)
names(Results_RF_plots_SR) <- c("Commu.AUC", "Max.Sorensen")

list_plots <- list(plots = vector(length = 911))
RF_Comparison_plots_SR <- matrix(ncol = 4, nrow = 911)
colnames(RF_Comparison_plots_SR) <- c("Plots", "CommunAUC", "MaxSorensen", "Model")
RF_Comparison_plots_SR[,1] <- seq(1,911)
RF_Comparison_plots_SR[,4] <- rep("Random.Forest.SR",911)

for (plots in 1:911) {
  
  if(mean(threeDMatrix[plots,,1]) == 0) {
    RF_Comparison_plots_SR[plots,2] <- NA }
  else {
    RF_Comparison_plots_SR[plots,2] <- auc(roc(threeDMatrix[plots,,10], as.factor(threeDMatrix[plots,,1]))) }
  
  RF_Comparison_plots_SR[plots,3] <- MaxSorensen(as.numeric(as.vector(append(threeDMatrix[plots,,1], threeDMatrix[plots,,10]))))
  
}

#########################################################################################################
####################################### GLM SPECIES RICHNESS ############################################
#########################################################################################################

Results_GLM5 <- vector(length = 4)
names(Results_GLM5) <- c("Mass.TSS", "Max.Kappa", "AUC", "Rsquared")

list_species <- list(species = vector(length = 100))
Comparison_GLM5 <- matrix(ncol = 6, nrow = 100)
colnames(Comparison_GLM5) <- c("Species", "TSS", "Kappa", "AUC", "Rsquared", "Model")
Comparison_GLM5[,1] <- seq(1,100)
Comparison_GLM5[,6] <- rep("GLM5",100)

Rsquared_Tjur_GLM5 <- vector(length = 10)

for (species in 1:100) {
  
  Comparison_GLM5[species,2] <- ecospat.max.tss(threeDMatrix[,species,11], threeDMatrix[,species,1])[[2]][1,2]
  Comparison_GLM5[species,3] <- ecospat.max.kappa(threeDMatrix[,species,11], threeDMatrix[,species,1])[[2]][1,2]
  Comparison_GLM5[species,4] <- auc(roc(threeDMatrix[,species,11], as.factor(threeDMatrix[,species,1])))
  
  for (trials in 1:10) {
    Rsquared_Tjur_GLM5[trials] <- PseudoR2(models_fitting[[10]][[trials]][[species]], which = "Tjur")
  }
  
  Comparison_GLM5[species,5] <- mean(Rsquared_Tjur_GLM5)
}

Results_GLM5_plots <- vector(length = 2)
names(Results_GLM5_plots) <- c("Commu.AUC", "Max.Sorensen")

list_plots <- list(plots = vector(length = 911))
Comparison_GLM5_plots <- matrix(ncol = 4, nrow = 911)
colnames(Comparison_GLM5_plots) <- c("Plots", "CommunAUC", "MaxSorensen", "Model")
Comparison_GLM5_plots[,1] <- seq(1,911)
Comparison_GLM5_plots[,4] <- rep("GLM5",911)

for (plots in 1:911) {
  
  if(mean(threeDMatrix[plots,,1]) == 0) {
    Comparison_GLM5_plots[plots,2] <- NA }
  else {
    Comparison_GLM5_plots[plots,2] <- auc(roc(threeDMatrix[plots,,11], as.factor(threeDMatrix[plots,,1]))) }
  
  Comparison_GLM5_plots[plots,3] <- MaxSorensen(as.numeric(as.vector(append(threeDMatrix[plots,,1], threeDMatrix[plots,,11]))))
  
}

#########################################################################################################
########################################## GRAPH SPECIES ADDITIONAL #####################################
#########################################################################################################

sp_additional_analyses <- rbind(RF_Comparison, RF_Comparison_SR, Comparison_GLM5)
sp_additional_analyses <- sp_additional_analyses[,-5]

glm4_tss <- lmer_analyses_df$TSS$Metrics[lmer_analyses_df$TSS$Model == "GLM4"]
glm4_kappa <- lmer_analyses_df$Kappa$Metrics[lmer_analyses_df$Kappa$Model == "GLM4"]
glm4_auc <- lmer_analyses_df$AUC$Metrics[lmer_analyses_df$AUC$Model == "GLM4"]

colnames(glm4_tss)[2] <- "TSS"
colnames(glm4_kappa)[2] <- "Kappa"
glm4_kappa <- glm4_kappa[,-3]
colnames(glm4_auc)[2] <- "AUC"
glm4_auc <- glm4_auc[,-3]

glm41 <- merge(glm4_tss, glm4_kappa, by = "Species")
glm4 <- merge(glm41, glm4_auc, by = "Species")

glm4$Species <- seq(1,100)

sp_additional_analyses <- rbind(sp_additional_analyses, glm4)
sp_additional_analyses$Model <- factor(sp_additional_analyses$Model, levels = c("GLM4", "GLM5", "Random.Forest.prev", "Random.Forest.SR"))

tss_ad <- ggplot(data = sp_additional_analyses, aes(x = Model, y = TSS, fill = Model)) +
  geom_boxplot() + labs(x = "Models", y = "Maximum TSS") +
  scale_fill_manual(values = c("GLM4" = "orange3", "Random.Forest.prev" = "plum3", "Random.Forest.SR" = "maroon4", "GLM5" = "firebrick3"),
                    breaks = c("GLM4", "Random.Forest.prev", "Random.Forest.SR", "GLM5"),
                    labels = c("Biotic GLM2", "Random Forest 2", "Random Forest SR", "Biotic GLM SR")) +
  scale_x_discrete(breaks = c("GLM4", "Random.Forest.prev", "Random.Forest.SR", "GLM5"),
                   labels = c("Biotic GLM2", "Random Forest 2", "Random Forest SR", "Biotic GLM SR"), name = NULL) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

kappa_ad <- ggplot(data = sp_additional_analyses, aes(x = Model, y = Kappa, fill = Model)) +
  geom_boxplot() + labs(x = "Models", y = "Maximum Kappa") +
  scale_fill_manual(values = c("GLM4" = "orange3", "Random.Forest.prev" = "plum3", "Random.Forest.SR" = "maroon4", "GLM5" = "firebrick3"),
                    breaks = c("GLM4", "GLM5", "Random.Forest.prev", "Random.Forest.SR"),
                    labels = c("Biotic GLM2", "Biotic GLM SR", "Random Forest 2", "Random Forest SR")) +
  scale_x_discrete(breaks = c("GLM4", "Random.Forest.prev", "Random.Forest.SR", "GLM5"),
                   labels = c("Biotic GLM2", "Random Forest 2", "Random Forest SR", "Biotic GLM SR"), name = NULL) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

auc_ad <- ggplot(data = sp_additional_analyses, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot() + labs(x = "Models", y = "AUC") +
  scale_fill_manual(values = c("GLM4" = "orange3", "Random.Forest.prev" = "plum3", "Random.Forest.SR" = "maroon4", "GLM5" = "firebrick3"),
                    breaks = c("GLM4", "Random.Forest.prev", "Random.Forest.SR", "GLM5"),
                    labels = c("Biotic GLM2", "Random Forest 2", "Random Forest SR", "Biotic GLM SR")) +
  scale_x_discrete(breaks = c("GLM4", "Random.Forest.prev", "Random.Forest.SR", "GLM5"),
                   labels = c("Biotic GLM2", "Random Forest 2", "Random Forest SR", "Biotic GLM SR"), name = NULL) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

grid.arrange(tss_ad, kappa_ad, auc_ad, nrow = 2, ncol = 2)

#########################################################################################################
########################################## TESTS SPECIES ADDITIONAL #####################################
#########################################################################################################
sp_additional_analyses$Species <- as.factor(sp_additional_analyses$Species)

lm_analyses_ad <- vector(mode = "list", length = 3)
Anova_analyses_ad <- vector(mode = "list", length = 3)
Metrics <- c("TSS", "Kappa", "AUC")

for (metrics in 1:3) {
  lm_analyses_ad[[metrics]] <- lm(paste(Metrics[metrics], " ~ Model + Species", sep = ""), data = sp_additional_analyses)
  Anova_analyses_ad[[metrics]] <- anova(lm_analyses_ad[[metrics]])
}

aov_analyses_ad <- vector(mode = "list", length = 3)
Tukey_analyses_ad <- vector(mode = "list", length = 3)
for (metrics in 1:3) {
  aov_analyses_ad[[metrics]] <- aov(as.formula(paste(Metrics[metrics], " ~ Model + Species", sep = "")), data = sp_additional_analyses)
  Tukey_analyses_ad[[metrics]] <- TukeyHSD(aov_analyses_ad[[metrics]])
}

Good_p_values <- vector(mode = "list", length = 3)
for (metrics in 1:3) {
  Good_p_values[[metrics]] <- Tukey_analyses_ad[[metrics]]$Model[,4][Tukey_analyses_ad[[metrics]]$Model[,4] < 0.05]
}


TSSP_values <- Tukey_analyses_ad[[1]]$Model[,4]
TSSP_values_corr <- TSSP_values*3

KappaP_values <- Tukey_analyses_ad[[2]]$Model[,4]
KappaP_values_corr <- KappaP_values*3

AUCP_values <- Tukey_analyses_ad[[3]]$Model[,4]
AUCP_values_corr <- AUCP_values*3

#########################################################################################################
############################################ GRAPH PLOTS ADDITIONAL #####################################
#########################################################################################################
plots_additional_analyses <- rbind(RF_Comparison_plots, RF_Comparison_plots_SR, Comparison_GLM5_plots)

glm3_AUCC <- lmer_analyses_df_plots$CommunAUC$Metrics[lmer_analyses_df_plots$CommunAUC$Model == "GLM3"]
glm3_MSI <- lmer_analyses_df_plots$Max_Sorensen$Metrics[lmer_analyses_df_plots$Max_Sorensen$Model == "GLM3"]

colnames(glm3_AUCC)[2] <- "CommunAUC"
colnames(glm3_MSI)[2] <- "MaxSorensen"
glm3_MSI <- glm3_MSI[,-3]

glm3 <- merge(glm3_AUCC, glm3_MSI, by = "Plots")

plots_additional_analyses <- rbind(plots_additional_analyses, glm3)
plots_additional_analyses$Model <- factor(plots_additional_analyses$Model, levels = c("GLM3", "GLM5", "Random.Forest.cor", "Random.Forest.SR"))

cAUC_ad <- ggplot(data = plots_additional_analyses, aes(x = Model, y = CommunAUC, fill = Model)) +
  geom_boxplot() + labs(x = "Models", y = "Community AUC") +
  scale_fill_manual(values = c("GLM3" = "dodgerblue4", "Random.Forest.cor" = "plum3", "Random.Forest.SR" = "maroon4", "GLM5" = "firebrick3"),
                    breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                    labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR")) +
  scale_x_discrete(breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                   labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR"), name = NULL) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

msi_ad <- ggplot(data = plots_additional_analyses, aes(x = Model, y = MaxSorensen, fill = Model)) +
  geom_boxplot() + labs(x = "Models", y = "Maximum Sorensen Index") +
  scale_fill_manual(values = c("GLM3" = "dodgerblue4", "Random.Forest.cor" = "plum3", "Random.Forest.SR" = "maroon4", "GLM5" = "firebrick3"),
                    breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                    labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR")) +
  scale_x_discrete(breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                   labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR"), name = NULL) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10))

grid.arrange(cAUC_ad, msi_ad, nrow = 1, ncol = 2)

#########################################################################################################
############################################ TESTS PLOTS ADDITIONAL #####################################
#########################################################################################################
plots_additional_analyses$Plots <- as.factor(plots_additional_analyses$Plots)

lm_analyses_plots_ad <- vector(mode = "list", length = 2)
Anova_analyses_plots_ad <- vector(mode = "list", length = 2)
Metrics <- c("CommunAUC", "MaxSorensen")

for (metrics in 1:2) {
  lm_analyses_plots_ad[[metrics]] <- lm(paste(Metrics[metrics], " ~ Model + Plots", sep = ""), data = plots_additional_analyses)
  Anova_analyses_plots_ad[[metrics]] <- anova(lm_analyses_plots_ad[[metrics]])
}

aov_analyses_plots_ad <- vector(mode = "list", length = 3)
Tukey_analyses_plots_ad <- vector(mode = "list", length = 3)
for (metrics in 1:2) {
  aov_analyses_plots_ad[[metrics]] <- aov(as.formula(paste(Metrics[metrics], " ~ Model + Plots", sep = "")), data = plots_additional_analyses)
  Tukey_analyses_plots_ad[[metrics]] <- TukeyHSD(aov_analyses_plots_ad[[metrics]])
}

Good_p_values <- vector(mode = "list", length = 3)
for (metrics in 1:3) {
  Good_p_values[[metrics]] <- Tukey_analyses_plots_ad[[metrics]]$Model[,4][Tukey_analyses_plots_ad[[metrics]]$Model[,4] < 0.05]
}


CommunAUC_values <- Tukey_analyses_plots_ad[[1]]$Model[,4]
CommunAUC_values_corr <- CommunAUC_values*2

MaxSoIn_values <- Tukey_analyses_plots_ad[[2]]$Model[,4]
MaxSoIn_values_corr <- MaxSoIn_values*2

#########################################################################################################
############################################ ELEV GRAPH ADDITIONAL ######################################
#########################################################################################################

head(plots_additional_analyses)
colnames(plots_additional_analyses)[1] <- "Plot"
Elevation_graph_ad <- merge(plots_additional_analyses, Plants_df_scaled[,c(1,9)], by = "Plot")
colnames(Elevation_graph_ad)[5] <- "DEM"

min_df <- vector(length = 911)
for (plots in 1:911) {
  min_df[plots] <- min(Elevation_graph_ad$MaxSorensen[Elevation_graph_ad$Plot == plots])
}

min_df <- cbind(min_df,seq(1,911))
colnames(min_df)[2] <- "Plot"

max_df <- vector(length = 911)
for (plots in 1:911) {
  max_df[plots] <- max(Elevation_graph_ad$MaxSorensen[Elevation_graph_ad$Plot == plots])
}

max_df <- cbind(max_df, seq(1,911))
colnames(max_df)[2] <- "Plot"

Elevation_graph_ad_mm <- merge(Elevation_graph_ad, min_df, by = "Plot")
Elevation_graph_ad_mm <- merge(Elevation_graph_ad_mm, max_df, by = "Plot")
colnames(Elevation_graph_ad_mm)[4] <- "Models"
colnames(Elevation_graph_ad_mm)[6] <- "Min"
colnames(Elevation_graph_ad_mm)[7] <- "Max"


elev_plot_ad <- ggplot(data = Elevation_graph_ad_mm, aes(x = DEM_vec, y = as.numeric(as.character(MaxSorensen)) , colour = Models)) +
  labs(x = "Elevation [m]", y = "Maximum Sorensen Index") + geom_point(size = 5) +
  geom_segment(aes(x = DEM_vec, xend = DEM_vec, y = Min, yend = Max), col = "black") +
  scale_color_manual(values = c("GLM3" = "dodgerblue4", "Random.Forest.cor" = "orchid4", "Random.Forest.SR" = "violetred3", "GLM5" = "firebrick3"),
                     breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                     labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR")) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text=element_text(size=20), axis.title=element_text(size=24), legend.position = "none")

#########################################################################################################
#################################### ELEVATION GRAPH SR CI ADD ##########################################
#########################################################################################################


pmod <- ggplot(data = Elevation_graph_ad, aes(x = DEM_vec, y = as.numeric(as.character(MaxSorensen)), colour = Model, fill = Model)) +
  labs(x = "Elevation [m]", y = "Max Sorensen Index") +
  geom_smooth(method = "loess", se = TRUE) + ylim(c(0,1.2)) +
  scale_color_manual(values = c("GLM3" = "dodgerblue4", "Random.Forest.cor" = "orchid4", "Random.Forest.SR" = "violetred3", "GLM5" = "firebrick3"),
                     breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                     labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR")) +
  scale_fill_manual(values = c("GLM3" = "dodgerblue4", "Random.Forest.cor" = "orchid4", "Random.Forest.SR" = "violetred3", "GLM5" = "firebrick3"),
                    breaks = c("GLM3", "Random.Forest.cor", "Random.Forest.SR", "GLM5"),
                    labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR"), guide = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=18),legend.position = "none")

pSR <-  ggplot(data = Elevation_graph_gathered_mm, aes(x = Elevation, y = sp_richness, colour = sp_richness)) +
  geom_smooth(data = Plant_df_scaled_ordered, aes(x = DEM_vec, y = sp_richness),
              method = "loess", col = "black", span = 0.6, se = TRUE) + 
  labs(x = "Elevation [m]", y = "Species Richness") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), 
        axis.text=element_text(size=16), axis.title=element_text(size=18))


g1 <- ggplot_gtable(ggplot_build(pmod))
g2 <- ggplot_gtable(ggplot_build(pSR))

pp <- c(subset(g1$layout, name == "panel", se = t:r))

g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]

ax$widths <- rev(ax$widths)

g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

ggsave("plot.pdf", g, width = 10, height = 6)

#########################################################################################################
######################################## ELEVATION GRAPH TESTS ADS ######################################
#########################################################################################################

anova(lm(MaxSorensen~Model+DEM_vec, data = Elevation_graph_ad))

#########################################################################################################
##############################         SPECIES RICHNESS Correlation         #############################
#########################################################################################################

Pred_Sp_R <- as.data.frame(cbind(Plants_df_scaled$Plot,Plants_df_scaled$DEM_vec, vector(length = 911), vector(length = 911), vector(length = 911), vector(length = 911), vector(length = 911)))
colnames(Pred_Sp_R) <- c("Plots", "Elevation", "Obs_SpR", "GLM3_SpR", "GLM5_SpR", "RF_cor_SpR", "RF_SR_SpR")
Pred_Sp_R$Plots <- as.factor(Pred_Sp_R$Plots)

Pred_Sp_R$Obs_SpR <- apply(threeDMatrix[,,1], MARGIN = 1, FUN = sum)
Pred_Sp_R$GLM3_SpR <- apply(threeDMatrix[,,6], MARGIN = 1, FUN = sum)
Pred_Sp_R$GLM5_SpR <- apply(threeDMatrix[,,8], MARGIN = 1, FUN = sum)
Pred_Sp_R$RF_cor_SpR <- apply(threeDMatrix[,,10], MARGIN = 1, FUN = sum)
Pred_Sp_R$RF_SR_SpR <- apply(threeDMatrix[,,11], MARGIN = 1, FUN = sum)

Pred_Sp_R_gathered <- gather(Pred_Sp_R, key = "Model", value = "Species_Richness", Obs_SpR:RF_SR_SpR)
Pred_Sp_R_gathered_2 <- gather(Pred_Sp_R, key = "Model", value = "Species_Richness", GLM3_SpR:RF_SR_SpR)
Pred_Sp_R_gathered_2$Model <- factor(Pred_Sp_R_gathered_2$Model, levels = c("GLM3_SpR", "GLM5_SpR", "RF_cor_SpR", "RF_SR_SpR"))

colnames(Pred_Sp_R_gathered_2)[4] <- "Models"

ggplot(Pred_Sp_R_gathered_2, aes(x = Obs_SpR, y = as.numeric(as.character(Species_Richness)), colour = Models, fill = Models)) + 
  geom_point() + labs(x = "Observed Species Richness", y = "Predicted Species Richness") +
  geom_abline(slope = 1, intercept = 0, col = "black") + xlim(c(0,50)) + ylim(c(0,50)) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_color_manual(values = c("GLM3_SpR" = "dodgerblue4", "RF_cor_SpR" = "plum3", "RF_SR_SpR" = "maroon4", "GLM5_SpR" = "firebrick3"),
                     breaks = c("GLM3_SpR", "RF_cor_SpR", "RF_SR_SpR", "GLM5_SpR"),
                     labels = c("Biotic GLM1", "Random Forest 1", "Random Forest SR", "Biotic GLM SR")) +
  scale_fill_manual(values = c("GLM3_SpR" = "dodgerblue4", "RF_cor_SpR" = "plum3", "RF_SR_SpR" = "maroon4", "GLM5_SpR" = "firebrick3"),
                    breaks = c("GLM3_SpR", "RF_cor_SpR", "RF_SR_SpR", "GLM5_SpR"), name = NULL, guide = FALSE) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 16), axis.text = element_text(size = 16), legend.position = "none")

#########################################################################################################
######################################### HMSC all data to choose 3 sp ##################################
#########################################################################################################

eco_var_tot <- Plants_df_scaled[,4:8]
secmat <- matrix(ncol = 5, nrow = nrow(Plants_df_scaled))
secmat <- as.data.frame(secmat)
secmat <- sapply(eco_var_tot, function(x) x^2)
colnames(secmat) <- c("GS_rad_sec", "GS_p_sec", "GS_t_sec", "pH_sec", "slope_sec")
var <- cbind(eco_var_tot, secmat)
# Creating a matrix for random factors
Random_factor <- matrix(ncol = 2, nrow = 911)
rownames(Random_factor) <- Plants_df_scaled$Plot
colnames(Random_factor) <- c("Random_1", "Random_elev")
# Filling the Random_factor$Random_1 column only with 1 for each plots (as they are all in the same
# area and all sampled the same way)
Random_factor[,1] <- rep(1,911)
#Filling the Random_factor$Random_elev according to the elvation class (Low or high)
for(rows in 1:nrow(Plants_df_scaled)){
  if(Plants_df_scaled$DEM_vec[rows] < median(Plants_df_scaled$DEM_vec))
    Random_factor[rows,2] <- 2
  else
    Random_factor[rows,2] <- 3
}
# Changing it into a dataframe and setting columns as factor
Random_factor <- as.data.frame(Random_factor)
Random_factor$Random_1 <- as.factor(Random_factor$Random_1)
Random_factor$Random_elev <- as.factor(Random_factor$Random_elev)
Test_all <- as.HMSCdata(Y = as.matrix(Plants_df_scaled[,10:109]), X = as.matrix(var), Random = Random_factor)


mod_test <- hmsc(Test_all, family = "probit", niter = 10000, nburn = 5000, thin = 50)

cor_Mat_Test <- corRandomEff(mod_test, cor = TRUE)

dimnames(cor_Mat_Test)[[1]] <- seq(1,100)
dimnames(cor_Mat_Test)[[2]] <- seq(1,100)

averageCor_Mat_Test <- apply(cor_Mat_Test[, , , 1], 1:2, mean)
corrplot(averageCor_Mat_Test, method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "lower")

species_corr_sum <- apply(abs(averageCor_Mat_Test), 1, sum)
species_corr_max <- order(species_corr_sum, decreasing = TRUE)
head(species_corr_max)

#########################################################################################################
############################################ Sign of INTERACTION ! ######################################
#########################################################################################################
All_int_neg <- vector(length = 100)
All_int_pos <- vector(length = 100)
# ALL SPECIES TO COMPARE
for (species in 1:100) {
  All_int_neg[species] <- sum(averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] < 0])
  All_int_pos[species] <- sum(averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] > 0])
}

mean(All_int_neg) # -5.10
mean(All_int_pos) # 6.12

# CORRELATION SPECIES
for (species in species_corr_max[1:3]) {
  Neg_int <- averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] < 0]
  Pos_int <- averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] > 0]
  print(paste("Species n° ", species, ", Negative interactions score = ", sum(Neg_int), ", Positive interactions score = ", sum(Pos_int), sep = ""))
}

# PREVALENT SPECIES
for (species in c(39, 94, 7)) {
  Neg_int <- averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] < 0]
  Pos_int <- averageCor_Mat_Test[,species][averageCor_Mat_Test[,species] > 0]
  print(paste("Species n° ", species, ", Negative interactions score = ", sum(Neg_int), ", Positive interactions score = ", sum(Pos_int), sep = ""))
}
