library(tidyverse)
library(ggplot2)

# Question 1
mu <- 0.15
V_Bar <- 0.04
delta_v <- 0.2
kappa_v <- 1
S0 <- 0.5
V0 <- V_Bar
Time_T <- 1
Step_Size <- 1000
N_Simulations <- 10000
dt <- Time_T/Step_Size
Corr <- c(-0.5, 0, 0.5)

Sim_SV <- function(rho, delta_v, seed = 42){
  set.seed(seed)
  
  S <- matrix(NA_real_, nrow = N_Simulations, ncol = Step_Size +1 )
  V <- matrix(NA_real_, nrow = N_Simulations, ncol = Step_Size +1 )
  
  S[,1] <- S0
  V[,1]<- V0
  
  for (t in 1:Step_Size) {
    z1 <- rnorm(N_Simulations)
    z2 <- rho* z1 + sqrt(1-rho^2) * rnorm(N_Simulations)
    
    V_T <- pmax(V[,t], 0)
    
    V[, t+1] <- pmax(
      V[, t] + kappa_v * (V_Bar - V[, t]) * dt +
        delta_v * sqrt(V_T) * sqrt(dt) * z2,
      0
    )

    S[, t+1] <- S[, t] * exp((mu - 0.5 * V_T) * dt + sqrt(V_T * dt) * z1)
    
  }
  
  tibble(
    Final_S = S[, Step_Size + 1],
    rho = factor(rho),
    delta_v = delta_v
  )
}

set.seed(123)

RhoTest <- map_dfr(Corr, ~Sim_SV(rho = .x, delta_v = 0.2))

RhoTest |>
  ggplot(aes(x = Final_S, fill = rho))+
  geom_density(alpha = 0.5) +
  scale_x_log10()+
  labs(
    title = "Log-Scale Density of S at t=1",
    x = "Log Scale",
    y = "Density"
  )

Delta_Array <- c(0.5, 1, 1.5)
New_RhoTest <- map_dfr(Delta_Array, ~Sim_SV(rho = 0, delta_v = .x))

New_RhoTest |>
  ggplot(aes(x = Final_S, fill = factor(delta_v)))+
  geom_density(alpha = 0.5)+
  labs(
    title = " Effect of δv on Final Stock Price",
    x = "Stock Price",
    y = "Density",
    fill = "δv"
  )

New_RhoTest |>
  ggplot(aes(x = Final_S, fill = factor(delta_v)))+
  geom_density( alpha = 0.5) +
  scale_x_log10()+
  labs(
    title = "Effect of δv ",
    x = "Log Scale",
    y = "Density",
    fill = "δv"
  )

lognorm_dist <- RhoTest |>
  group_by(rho) |>
  summarise(
    mu_hat = mean(log(Final_S)),
    sigma_hat = sd(log(Final_S)),
    .groups = "drop"
  )

lognorm_dist

Difference <- RhoTest |>
  group_by(rho) |>
  summarise(x = list(seq(min(Final_S), quantile(Final_S, 0.999), length.out = 500)))|>
  unnest(x) |>
  left_join(lognorm_dist, by = "rho") |>
  mutate(dln = dlnorm(x, meanlog = mu_hat, sdlog = sigma_hat))
  
ggplot(RhoTest, aes(x = Final_S, color = rho)) +
  geom_density(aes(fill = rho), alpha = 0.3)+
  geom_line(data = Difference, aes(x=x, y = dln, color = rho),
            linewidth = 1) +
  scale_x_log10()+
  labs(
    title = "Empirical Vs Log-Normal Densities",
    x="Log Scale",
    y = "Density"
  )
  
#Question 2
# a) GBM Discretization
S0 <- 100
r <- 0.05
sigma <- 0.2
Time_T <- 1
N <- 50
M <- 100000
dt <- Time_T/N
Simm <- (1:N)*dt

set.seed(1)

GBM <- function(M,N,dt,S0, r, sigma){
  z <- matrix(rnorm(M*N), nrow = M, ncol = N)
  incr <- (r - 0.5* sigma^2) * dt + sigma * sqrt(dt) * z
  log_s <- cbind(log(S0), t(apply(incr, 1,cumsum)))
  S <- exp(log_s)
  S
}

paths <- GBM(M,N,dt,S0,r,sigma)
S_grid <- paths[, -1, drop = FALSE]
S_T <- S_grid[, N]
S_Bar <- rowMeans(S_grid)

payoff_x <- pmax(S_T - S_Bar, 0)
discrete_x <- exp(-r * Time_T) * payoff_x

Cmc <- mean(discrete_x)
SE_Cmc <- sd(discrete_x)/sqrt(M) #Standard Error
CI_Cmc <- Cmc + c(-1,1) * 1.96 * SE_Cmc # Confidence Interval

list(
  CMC = Cmc,
  CI95 = CI_Cmc
)

# b) Control Variate
Asian_Call <- function(S0, r, sigma, Time_T, N){
  ti <- (1:N)*(Time_T/N)
  Summ_t <- sum(ti)
  Summ_min <- sum(pmin(rep(ti, each =N), rep(ti, times = N)))
  
  M_G <- log(S0) + (r - 0.5 * sigma^2) * (Summ_t/N)
  V_G <- (sigma^2 / N^2) * Summ_min
  
 Var_ST <- sigma^2 * Time_T
 Cov <- (sigma^2/N) * Summ_t
 Var_ex <- Var_ST + V_G -2 * Cov
 sigma_ex <- sqrt(Var_ex/Time_T)
 
 S1 <- S0
 S2 <- exp(M_G + 0.5 * V_G - r * Time_T)
 
 d1 <- (log(S1/S2)+0.5 * sigma_ex^2 * Time_T)/(sigma_ex* sqrt(Time_T))
 d2 <- d1 - sigma_ex * sqrt(Time_T)
 
 mu_Y <- S1 * pnorm(d1) - S2 * pnorm(d2)
 mu_Y
}

mu_Y <- Asian_Call(S0, r, sigma, Time_T, N)
Geo_S <- exp(rowMeans(log(S_grid)))
payoff_y <- pmax(S_T - Geo_S, 0)
discrete_y <- exp(-r * Time_T) * payoff_y

theta_hat <-cov(discrete_x, discrete_y) / var(discrete_y)
z <- discrete_x+ theta_hat * (mu_Y - discrete_y)

Ccv <- mean(z)
SE_Ccv <- sd(z)/sqrt(M)
CI_Ccv <- Ccv + c(-1,1) * 1.96 * SE_Ccv

list(
  CCV = Ccv,
  CI95_Ccv = CI_Ccv,
  theta_hat = theta_hat,
  mu_Y = mu_Y
)

# c)
Var_MC <- var(discrete_x)
Var_CV <- var(z)
VR <- 100*(1 - Var_CV/Var_MC)

tibble(
  estimator = c("Plain MC", "With Control Variate"),
  price = c(Cmc, Ccv),
  SE = c(SE_Cmc, SE_Ccv),
  CI_Low = c(CI_Cmc[1], CI_Ccv[1]),
  CI_High = c(CI_Cmc[2], CI_Ccv[2])
) |>
  print()

cat(sprintf("Variance Reduction: %.1f%%\n", VR)) # Variance Reduction: 99%, 0.1% of the variance is uncorrelated to the plain MC

# d) CI 
CV_Price <- function(M,N,sigma, seed=1){
  set.seed(seed)
  dt<- Time_T/N
  S_Paths <- GBM(M,N,dt,S0,r,sigma)
  S_grid <- S_Paths[,-1,drop = FALSE]
  S_T <- S_grid[,N]
  S_Bar <- rowMeans(S_grid)
  X<- exp(-r * Time_T) * pmax(S_T -S_Bar,0)
  
  G <- exp(rowMeans(log(S_grid)))
  Y <- exp(-r *Time_T) * pmax(S_T - G, 0)
  MuY <- Asian_Call(S0, r, sigma, Time_T, N)
  
  theta <- cov(X,Y)/var(Y)
  z <- X + theta * (MuY - Y)
  
  c(MC = mean(X), SE_MC = sd(X)/sqrt(M),
  CV = mean(z), SE_Cv = sd(z)/sqrt(M)
  )
}

Ms <- round(exp(seq(log(1e3), log(1e5), length.out = 6)))
res_m <- map_dfr(Ms, \(m) {
  out <- CV_Price(M = m, N = N, sigma = sigma, seed = 123)
  tibble(
    M = m,
    mc_low = out["MC"] - 1.96 * out["SE_MC"],
    mc_high = out["MC"] + 1.96 * out["SE_MC"],
    cv_low = out["CV"] - 1.96 * out["SE_Cv"],
    cv_high = out["CV"] + 1.96 * out["SE_Cv"]
  )
})

res_m |>
  mutate(
    width_mc = mc_high - mc_low,
    width_cv = cv_high - cv_low
  ) |>
  pivot_longer(starts_with("width_"), names_to = "method", values_to = "width") |>
  mutate(method = recode(method, width_mc = "Plain MC", width_cv = "Control variate")) |>
  ggplot(aes(x = M, y = width, color = method)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  labs(title = "95% CI width vs simulations M",
       x = "log scale", y = "CI width") +
  theme_minimal()



