agalinis.height.long$week <- as.numeric(gsub("t","", agalinis.height.long$week))

fit_logistic_model2 <- function(agalinis.height.long) {
  time <- agalinis.height.long$week
  height <- agalinis.height.long$height
  
  logistic_growth <- function(time,r,maxH) {
    N0 <- height[1]
    N <- N0 * maxH / (N0 + (maxH - N0) * exp(-r * time))
    return(N)
  }
  
  #MLE 
  model <- optim(par = c(r = 0.05, maxH = 40),
                 fn = function(par) sum((height - logistic_growth(time,par[1],par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  maxH_hat <- model$par[2]
  
  #predictions
  predicted_height <- logistic_growth(time,r_hat,maxH_hat)
  
  return(data.frame(species = unique(agalinis.height.long$species),
                    type = unique(agalinis.height.long$type),
                    treatment = unique(agalinis.height.long$treatment),
                    week = time,
                    observed_height = height,
                    predicted_height = predicted_height,
                    r = r_hat,
                    maxH = maxH_hat,
                    N0 = height[1]))
}


results_height_agalinis <- agalinis.height.long %>%
  filter(!is.na(height)) %>%
  group_by(treatment,species,type,replicate_id) %>%
  do(fit_logistic_model2(.))

results_height_agalinis %>%
  mutate(unique_ID = paste(species,treatment,type,replicate_id,sep="_")) %>%
  ggplot(aes(x = week, y = predicted_height)) +
  geom_line() +
  facet_wrap(~unique_ID,scales="free") +
  geom_point(aes(x=week,y=observed_height),size=2,alpha=0.3) +
  labs(x="Week",y="Height") +
  theme_classic()

ggplot(results_height_agalinis,aes(x=type,y=r)) +
  geom_boxplot() +
  theme_classic()

# try with just HESU
hetero.height.long$week <- as.numeric(gsub("t","", hetero.height.long$week))

fit_logistic_model3 <- function(hetero.height.long) {
  time <- hetero.height.long$week
  height <- hetero.height.long$height
  
  logistic_growth <- function(time,r,maxH) {
    N0 <- height[1]
    N <- N0 * maxH / (N0 + (maxH - N0) * exp(-r * time))
    return(N)
  }
  
  #MLE 
  model <- optim(par = c(r = 0.001, maxH = 14),
                 fn = function(par) sum((height - logistic_growth(time,par[1],par[2]))^2),
                 method = "BFGS",
                 control = list(maxit = 200))
  
  r_hat <- model$par[1]
  maxH_hat <- model$par[2]
  
  #predictions
  predicted_height <- logistic_growth(time,r_hat,maxH_hat)
  
  return(data.frame(species = unique(hetero.height.long$species),
                    type = unique(hetero.height.long$type),
                    treatment = unique(hetero.height.long$treatment),
                    week = time,
                    observed_height = height,
                    predicted_height = predicted_height,
                    r = r_hat,
                    maxH = maxH_hat,
                    N0 = height[1]))
}

# fit whole dataset to get potential values for optimizer
# hetero.height.long %>%
#   filter(!is.na(height)) %>%
#   do(fit_logistic_model3(.)) %>%
#   summarize(mean_r = mean(r),
#             max_maxH = max(maxH))

results_height_hetero <- hetero.height.long %>%
  filter(!is.na(height)) %>%
  group_by(treatment,species,type,replicate_id) %>%
  do(fit_logistic_model3(.))

results_height_hetero %>%
  mutate(unique_ID = paste(species,treatment,type,replicate_id,sep="_")) %>%
  ggplot(aes(x = week, y = predicted_height)) +
  geom_line() +
  facet_wrap(~unique_ID,scales="free") +
  geom_point(aes(x=week,y=observed_height),size=2,alpha=0.3) +
  labs(x="Week",y="Height") +
  theme_classic()


# calculate relative growth rate
# add lag column for each individual
test <- height.long %>% 
  group_by() %>%
  mutate(height_t2 = dplyr::lead(height))

