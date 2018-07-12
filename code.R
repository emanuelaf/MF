rm(list=ls())

source("functions2.R")
set.seed(1)

data <- generate_pop(mu_domains = c(30, 30, 30, 30, 30, 30, 30), 
                     sd_domains = c(0.2, 0.2, 10, 0.2, 10, 10, 10),
                     frame_cost = c(6, 3, 1))

# print frame specific statistics
data %>% group_by(domain) %>% summarise(sd(Y))
data %>% group_by(domain) %>% summarise(mean(Y))
data %>% filter(A == 1) %>% summarise(mean_A = mean(Y), sd_A = sd(Y), cv_A = sd(Y)/mean(Y))
data %>% filter(B == 1) %>% summarise(mean_B =mean(Y), sd_B = sd(Y), cv_B = sd(Y)/mean(Y))
data %>% filter(C == 1) %>% summarise(mean_C = mean(Y), sd_C = sd(Y), cv_C = sd(Y)/mean(Y))

n <- 600

N <- nrow(data)

N_a <- data %>% filter(domain == "a") %>% nrow()
N_b_ab <- data %>% filter(domain == "b" | domain == "ab") %>% nrow()
N_C <- data %>% filter(domain == "c" | domain == "ac" | domain == "bc" | domain == "abc") %>% nrow()

N_A <- data %>% filter(domain == "a" | domain == "ab" | domain == "ac" | domain == "abc") %>% nrow()
N_B <- data %>% filter(domain == "b" | domain == "ab" | domain == "bc" | domain == "abc") %>% nrow()

# sample

n.sim <- 1000

res_mf_m <- numeric(n.sim)
res_mf_ka <- numeric(n.sim)
res_strat_reduced <- numeric(n.sim)
res_strat_new_size_1 <- numeric(n.sim)
res_strat_new_size_2 <- numeric(n.sim)
cost_s_scr_reduced_size <- numeric(n.sim)
cost_s_scr_new_size_1 <- numeric(n.sim)
cost_s_scr_new_size_2 <- numeric(n.sim)
cost_s_mf <- numeric(n.sim)
length_s_A_excluded <- numeric(n.sim)
length_s_B_excluded <- numeric(n.sim)


for (i in 1:n.sim)
{
# sample size allocation
  sigma2_alpha <- sigma_alpha_values(data = data, N_A = N_A, N_B = N_B, N_C = N_C)
  n <- n_size(type = "optimal", N = N, N_A = N_A, N_B = N_B, N_C = N_C, C = 1000, 
              c_0 = 100, sigma2_alpha = sigma2_alpha, n = 900)

  s_scr_reduced_size <- sample_scr_reduced_size(data, 
                                              n_A = n$n_A, 
                                              n_B = n$n_B,
                                              n_C = n$n_C)
  
  s_scr_new_size_1 <- sample_scr_new_size_1(data, 
                                            s_scr_reduced_size, 
                                            n_A = n$n_A, 
                                            n_B = n$n_B,
                                            n_C = n$n_C)
  
  s_scr_new_size_2 <- sample_scr_new_size_2(data,
                                            n_A = n$n_A, 
                                            n_B = n$n_B,
                                            n_C = n$n_C)
  
  s_mf <- s_scr_reduced_size
  s_mf$s_A_final <- s_scr_reduced_size$s_A_init
  s_mf$s_B_final <- s_scr_reduced_size$s_B_init
  s_mf$s_C_final <- s_scr_reduced_size$s_C_init

# number of rejections
  length_s_A_excluded[i] <- nrow(s_scr_reduced_size$s_A_excluded)
  length_s_B_excluded[i] <- nrow(s_scr_reduced_size$s_B_excluded)

# costs
# ampiezza campionaria di fatto dimezzata con metodo reduced_size
  cost_s_scr_reduced_size[i] <- cost(s_scr_reduced_size)
  cost_s_scr_new_size_1[i] <- cost(s_scr_new_size_1)
  cost_s_scr_new_size_2[i] <- cost(s_scr_new_size_2)
  cost_s_mf[i] <- cost(s_mf)

# estimators
  res_mf_m[i] <- est_mf_multiplicity(s_mf, N_A = N_A, N_B = N_B, N_C=N_C)
  res_mf_ka[i] <- est_mf_ka(s_mf, N_A = N_A, N_B = N_B, N_C = N_C)
  res_strat_reduced[i] <- est_strat(s_scr_reduced_size, N_a = N_a, N_b_ab = N_b_ab, N_C = N_C)
  res_strat_new_size_1[i] <- est_strat(s_scr_new_size_1, N_a = N_a, N_b_ab = N_b_ab, N_C = N_C)
  res_strat_new_size_2[i] <- est_strat(s_scr_new_size_2, N_a = N_a, N_b_ab = N_b_ab, N_C = N_C)
}

# risultati in boxplot
boxplot(data.frame(res_mf_m, res_mf_ka, res_strat_reduced, res_strat_new_size_1, res_strat_new_size_2))

data.frame(res_mf_m, res_strat_reduced, res_strat_new_size_1, res_strat_new_size_2) %>% 
  summarise(sd(res_mf_m), sd(res_mf_ka), sd(res_strat_reduced), sd(res_strat_new_size_1), sd(res_strat_new_size_2))

data.frame(cost_s_scr_reduced_size,
           cost_s_scr_new_size_1,
           cost_s_scr_new_size_2,
           cost_s_mf) %>% 
  summarise(mean(cost_s_mf), mean(cost_s_scr_reduced_size), mean(cost_s_scr_new_size_1), mean(cost_s_scr_new_size_2))

data.frame(length_s_A_excluded,
           length_s_B_excluded) %>% 
  summarise(min(length_s_A_excluded), 
            min(length_s_B_excluded),
            median(length_s_A_excluded), 
            median(length_s_B_excluded),
            mean(length_s_A_excluded), 
            mean(length_s_B_excluded),
            max(length_s_A_excluded), 
            max(length_s_B_excluded))

(var_mf_multiplicity(data, n_A = 100, n_B = 200, n_C = 300) - var(res_mf_m))/var_mf_multiplicity(data, n_A = 100, n_B = 200, n_C = 300)*100

