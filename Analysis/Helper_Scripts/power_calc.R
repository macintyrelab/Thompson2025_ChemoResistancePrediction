power_calc <- function(alpha=0.05, power=0.8, P=1, Rsquared=0, pred_prop, hr){
  # This function performs a power calculation specific to cox models, giving
  # as output the number of patients required to be powered.
  # Inputs: 
  #         alpha - the risk of rejecting a true null hypothesis, normally 0.05
  #         power - the likelihood of detecting an effect if it exists, normally 0.2
  #         P - the proportion of patients that did not have an event (e.g. death)
  #             before the end of the clinical history
  #         Rsquared - the Rsquared correlation between the predictor variable
  #                    and the covariates
  #         pred_prop - the proportion of patients who were predicted positive
  #         hr - the expected or observed hazard ratio
  # Outputs:
  #         d - the number of patients required in a powered cohort
  zscore_alpha <- qnorm(p=1-alpha, mean=0, sd=1) # Assume 1-tailed
  zscore_power <- qnorm(p=power, mean=0, sd=1)
  numerator <- (zscore_alpha + zscore_power)^2
  variance <- var(c(rep(1, pred_prop*10000),
                  rep(0, (1-pred_prop)*10000))) # Variance tends to a limit as n 
           # increases, so we approximate the limit here by multiplying by 10000
  Bsquared <- log(hr)^2
  denominator <- P * (1 - Rsquared) * variance * Bsquared
  d <- numerator / denominator
  d
}

power_calc_inv <- function(alpha=0.05, power=0.8, P=1, Rsquared=0, pred_prop, d){
  # This function functions similarly to power_calc, but instead for a given 
  # sample size d, it calculates the HR that the cohort is powered to detect.
  # Inputs: 
  #         alpha - the risk of rejecting a true null hypothesis, normally 0.05
  #         power - the likelihood of detecting an effect if it exists, normally 0.2
  #         P - the proportion of patients that did not have an event (e.g. death)
  #             before the end of the clinical history
  #         Rsquared - the Rsquared correlation between the predictor variable
  #                    and the covariates
  #         pred_prop - the proportion of patients who were predicted positive
  #         d - the number of patients in the cohort
  # Outputs:
  #         hr_upper - the upper HR that the cohort is powered to detect
  #         hr_lower - the lower HR that the cohort is powered to detect
  zscore_alpha <- qnorm(p=1-alpha, mean=0, sd=1) # Assume 1-tailed
  zscore_power <- qnorm(p=power, mean=0, sd=1)
  numerator <- zscore_alpha + zscore_power
  variance <- var(c(rep(1, pred_prop*10000),
                    rep(0, (1-pred_prop)*10000))) # Variance tends to a limit as n 
  # increases, so we approximate the limit here by multiplying by 10000
  denominator <- sqrt(P * (1 - Rsquared) * variance * d)
  B <- numerator / denominator
  hr <- sort(c(exp(B), 1/exp(B)))
  names(hr) <- c('lower', 'upper')
  hr
}