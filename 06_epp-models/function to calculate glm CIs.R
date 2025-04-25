
################################################################################
# Function for calculating 95% confidence intervals for GLMs
# Heather Kenny-Duddela
# Nov 25, 2024
################################################################################


calculate.glm.ci <- function(model, newdata) {
  # make a dataframe to hold the new data and results
  output.vals <- newdata
  
  # get predictions on the data scale
  output.vals$fit <- predict.glm(
    model,
    newdata, type="response")
  
  # get predictions on the link scale (scale of the model)
  link <- predict.glm(model, newdata=newdata,
                                   type="link", se.fit=T)
  
  # save the fitted values
  output.vals$fit_link <- link$fit
  output.vals$se_link <- link$se.fit
  # get inverse-link function
  fam.model <- family(model)
  ilink.model <- fam.model$linkinv
  # calculate CIs as back-transformed fit_link +/- 2*se_link
  output.vals$upper <- ilink.model(output.vals$fit_link + 
                                              2*output.vals$se_link)
  output.vals$lower <- ilink.model(output.vals$fit_link - 
                                              2*output.vals$se_link)
  return(output.vals)
}


