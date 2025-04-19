library(tidyverse)
library(survival)
library(ggplot2)
valung <- read.csv("C:/Users/86181/Downloads/valung.csv")

valung <- valung %>%
  mutate(
    t = ifelse(t <= 0, 0.1, t),     
    log_time = log(t),
    dead = ifelse(dead == "dead", 1, 0), 
    therapy = factor(therapy),
    cell = factor(cell),
    prior = factor(prior, levels = c("no", "yes"))
  )


summary(valung)
table(valung$therapy)
table(valung$cell)

ggplot(valung, aes(x = t, fill = as.factor(dead))) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.7) +
  labs(title = "Survival Time Distribution", x = "Time (days)", fill = "Status")

surv_obj <- Surv(time = valung$t, event = valung$dead)
fit_km <- survfit(surv_obj ~ therapy, data = valung)

valung <- valung %>%
  mutate(
    cell_squamous = ifelse(cell == "Squamous", 1, 0),
    cell_small = ifelse(cell == "Small", 1, 0),
    cell_adeno = ifelse(cell == "Adeno", 1, 0),
    cell_large = ifelse(cell == "Large", 1, 0),
    therapy_bin = ifelse(therapy == "test", 1, 0),
    prior_bin = ifelse(prior == "yes", 1, 0)
  )

X <- valung %>%
  select(therapy_bin, age, kps, diagtime, prior_bin,
         cell_small, cell_adeno, cell_large) %>%
  as.matrix()

X_scaled <- scale(X)

stan_data <- list(
  N = nrow(X_scaled),
  D = ncol(X_scaled),
  X = X_scaled,
  y = valung$log_time,
  status = valung$dead
)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_fit <- stan(
  file = "C:\\Users\\86181\\Desktop\\STAT 447\\STAT447C_project\\gp_survival_model.stan",  
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 1234,
  init = 0,
  control = list(adapt_delta = 0.95)
)
traceplot(stan_fit, pars = c("alpha", "sigma", "ell[1]", "ell[2]", "lp__"))

summary(stan_fit)$summary[,"Rhat"]

posterior_samples <- extract(stan_fit)

alpha_samples <- posterior_samples$alpha
hist(alpha_samples, main = "Posterior of alpha", xlab = "alpha")

f_hat <- colMeans(posterior_samples$f)
plot(valung$log_time, f_hat,
     xlab = "Observed log(time)",
     ylab = "Predicted log(time)",
     main = "Posterior Mean Predictions")
abline(0, 1, col = "red", lty = 2)

y_rep <- posterior_samples$y_rep

test_idx <- which(valung$therapy_bin == 1)
standard_idx <- which(valung$therapy_bin == 0)

exp_test <- exp(rowMeans(y_rep[, test_idx]))
exp_standard <- exp(rowMeans(y_rep[, standard_idx]))

mean_diff <- exp_test - exp_standard
mean(mean_diff)
quantile(mean_diff, c(0.025, 0.975))

ell_samples <- posterior_samples$ell

ell_summary <- apply(ell_samples, 2, function(x) c(mean = mean(x),
                                                   q025 = quantile(x, 0.025),
                                                   q975 = quantile(x, 0.975)))
ell_summary <- t(ell_summary)
colnames(ell_summary) <- c("Mean", "2.5%", "97.5%")


colnames(X_scaled)  
rownames(ell_summary) <- colnames(X_scaled)

ell_summary[order(ell_summary[, "Mean"]), ]

