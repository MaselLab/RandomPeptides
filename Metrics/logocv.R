# Leave one group out cross validation.

# This function assumes that your group is your random intercept.

logo_cv <- function(
  data.df,
  predictors, # List of all predictors as strings.
  dependent.variable, # String of the dependent variable.
  group.col,
  weight.col,
  intercept = T,
  verbose = T,
  output = "WMSE", # Currently can be set to WMSE. If blank, the whole DF with pred column will be returned.
  ...
){
  require(lme4)
  require(tidyverse)
  require(Hmisc)
  
  # Adding columns to data frame.
  validation.df <- tibble(
    "Y" = data.df[[dependent.variable]],
    "Weights" = data.df[[weight.col]],
    "Groups" = data.df[[group.col]],
    "Pred" = NA
  )
  
  for (i in 1:length(predictors)) {
    validation.df[[paste0("X", i)]] <- data.df[[predictors[i]]]
  }
  
  if (verbose == T){
    print("Data frame for analysis:")
    print(validation.df)
  }
  
  # Setting up the groups and variables for the LOGO CV.
  groups <- unique(validation.df$Groups)
  
  model.lmer <- NA
  
  pred.lmer <- NA
  
  # Writing the model formula.
  names.list <- names(validation.df)
  model.formula <-
    as.formula(
      paste("Y ~", paste(names.list[!names.list %in% c("Y", "Weights", "Groups", "Pred")], collapse = " + "), "+ (1|Groups)")
    )
  
  if (verbose == T){
    print("Formula for LOGOCV:")
    print(model.formula)
  }
  
  for (i in 1:length(groups)) {
    train.data <- validation.df %>% filter(Groups != groups[i])
    test.data <- validation.df %>% filter(Groups == groups[i])
    
    #print(train.data)
    #print(test.data)
    
    model.lmer <- lmer(
      formula = model.formula,
      data = train.data,
      weights = Weights
    )
    
    pred.lmer <- predict(model.lmer, newdata = test.data, ...)
    
    validation.df[validation.df$Groups == groups[i],]$Pred <- pred.lmer
  }
  
  if (output == "WMSE"){
    wmse.df <- validation.df %>% group_by(Groups) %>%
      summarise(
        W_group = sum(Weights),
        Pred_group = wtd.mean(Pred, weights = Weights),
        Y_group = wtd.mean(Y, weights = Weights)
      )
    wmse.df$Res_group <- wmse.df$Y_group - wmse.df$Pred_group
    wmse_numerator <-
      sum(wmse.df$W_group * (wmse.df$Res_group ^ 2))
    wmse_denominator <-
      sum(wmse.df$W_group)
    return(wmse_numerator / wmse_denominator)
  } else {
    return(validation.df)
  }
}