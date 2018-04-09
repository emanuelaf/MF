rm(list=ls())

source("functions2.R")

set.seed(1)
data <- generate_pop(shape_pars_frame = c(10, 5, 1), c(6, 4, 2))

n <- 30

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


for (i in 1:n.sim)
{
s_scr_reduced_size <- sample_scr_reduced_size(data, 
                                              n_A = proportional_n_screener(N=N, N_a = N_a, 
                                                                            N_b_ab = N_b_ab, 
                                                                            N_C = N_C, n = n)$n_a, 
                                              n_B = proportional_n_screener(N=N, N_a = N_a, 
                                                                            N_b_ab = N_b_ab, 
                                                                            N_C = N_C, n = n)$n_b_ab,
                                              n_C = proportional_n_screener(N=N, N_a = N_a, 
                                                                            N_b_ab = N_b_ab, 
                                                                            N_C = N_C, n = n)$n_C)

s_scr_new_size_1 <- sample_scr_new_size_1(data, 
                                          s_scr_reduced_size, 
                                          n_A = proportional_n_screener(N=N, N_a = N_a, 
                                                                        N_b_ab = N_b_ab, 
                                                                        N_C = N_C, n = n)$n_a, 
                                          n_B = proportional_n_screener(N=N, N_a = N_a, 
                                                                        N_b_ab = N_b_ab, 
                                                                        N_C = N_C, n = n)$n_b_ab,
                                          n_C = proportional_n_screener(N=N, N_a = N_a, 
                                                                        N_b_ab = N_b_ab, 
                                                                        N_C = N_C, n = n)$n_C)

s_scr_new_size_2 <- sample_scr_new_size_2(data,
                                          n_A = proportional_n_screener(N=N, N_a = N_a,                                     N_b_ab = N_b_ab, N_C = N_C, n = n)$n_a, 
                                          n_B = proportional_n_screener(N=N, N_a = N_a, 
                                                                        N_b_ab = N_b_ab, 
                                                                        N_C = N_C, n = n)$n_b_ab,
                                          n_C = proportional_n_screener(N=N, N_a = N_a, 
                                                                        N_b_ab = N_b_ab, 
                                                                        N_C = N_C, n = n)$n_C)

s_mf <- sample_mf(data,
                  n_A = proportional_n_mf(N=N, N_A = N_A, N_B = N_B, N_C = N_C, n = n)$n_A, 
                  n_B = proportional_n_mf(N=N, N_A = N_A, N_B = N_B, N_C = N_C, n = n)$n_B,
                  n_C = proportional_n_mf(N=N, N_A = N_A, N_B = N_B, N_C = N_C, n = n)$n_C)


# costs
# ampiezza campionaria di fatto dimezzata con metodo reduced_size
cost(s_scr_reduced_size)
cost(s_scr_new_size_1)
cost(s_scr_new_size_2)
cost(s_mf)


# estimators
res_mf_m[i] <- est_mf_multiplicity(s_mf)
res_mf_ka[i] <- est_mf_ka(s_mf, N_A = N_A, N_B = N_B, N_C=N_C)
res_strat_reduced[i] <- est_strat(s_scr_reduced_size)
res_strat_new_size_1[i] <- est_strat(s_scr_new_size_1)
res_strat_new_size_2[i] <- est_strat(s_scr_new_size_2)
}

par(mfrow = c(1, 5))
boxplot(res_mf_m)
boxplot(res_mf_ka)
boxplot(res_strat_reduced)
boxplot(res_strat_new_size_1)
boxplot(res_strat_new_size_2)

