require(dplyr)

### POPULATION

# generate pop with different variance differntials
# at the moment, N fixed and overlap fixed!!

generate_pop <- function(mu_domains, sd_domains, frame_cost) {
  N <- 70000
  
  # modify overlap here
  A <- c(rep(1, 40000), rep(0, 30000))
  B <- c(rep(0, 10000), rep(1, 10000), rep(0, 10000), rep(1, 30000), rep(0, 10000))
  C <- c(rep(0, 20000), rep(1, 30000), rep(0, 10000), rep(1, 10000))
  
  population <- data.frame(
    id = 1:N,
    A = A,
    B = B,
    C = C
    )

  population <- population %>% mutate(domain = ifelse(A == 1 & B == 0 & C == 0, "a",
                                          ifelse(A == 0 & B == 1 & C == 0, "b",
                                                 ifelse(A == 0 & B == 0 & C == 1, "c",
                                                        ifelse(A == 1 & B == 1 & C == 0, "ab",
                                                               ifelse(A == 1 & B == 0 & C == 1, "ac",
                                                                      ifelse(A == 0 & B == 1 & C == 1, "bc",
                                                                             ifelse(A == 1 & B == 1 & C == 1, "abc", NA))))
                                                 )))) %>%
    mutate(frame_cost_A = frame_cost[1], frame_cost_B = frame_cost[2], frame_cost_C = frame_cost[3]) %>%
    as_tibble()
  
  N_a <- length(population$domain[population$domain == "a"])
  N_b <- length(population$domain[population$domain == "b"])
  N_c <- length(population$domain[population$domain == "c"])
  N_ab <- length(population$domain[population$domain == "ab"])
  N_bc <- length(population$domain[population$domain == "bc"])
  N_ac <- length(population$domain[population$domain == "ac"])
  N_abc <- length(population$domain[population$domain == "abc"])

  
  population <- population %>% mutate(
    Y = ifelse(domain == "a", mu_domains[1] + rnorm(N_a, 0, sd_domains[1]),  
               ifelse(domain == "b", mu_domains[2] + rnorm(N_b, 0, sd_domains[2]),
               ifelse(domain == "c", mu_domains[3] + rnorm(N_c, 0, sd_domains[3]),
               ifelse(domain == "ab", mu_domains[4] + rnorm(N_ab, 0, sd_domains[4]),
               ifelse(domain == "bc", mu_domains[5] + rnorm(N_bc, 0, sd_domains[5]),
               ifelse(domain == "ac", mu_domains[6] + rnorm(N_ac, 0, sd_domains[6]),
               ifelse(domain == "abc", mu_domains[7] + rnorm(N_abc, 0, sd_domains[7]), 0)))))))
  )
  
  return(population)

}


### SCREENER DESIGN

# units not reallocated
sample_scr_reduced_size <- function(data, n_A, n_B, n_C) {
  s_A_init <- sample_n(data[data$A == 1 , ], n_A)
  s_A_final <- s_A_init %>% filter(domain == "a")
  s_A_excluded <- anti_join(s_A_init, s_A_final, by = c("id"))
  s_B_init <- sample_n(data[data$B == 1 , ], n_B)
  s_B_final <- s_B_init %>% filter(domain == "b" | domain == "ab")
  s_B_excluded <- anti_join(s_B_init, s_B_final, by = c("id"))
  s_C_init <- sample_n(data[data$C == 1 , ], n_C)
  s_C_final <- s_C_init 
  samples <- list(s_A_init = s_A_init, s_A_final = s_A_final, s_A_excluded = s_A_excluded,
                  s_B_init = s_B_init, s_B_final = s_B_final, s_B_excluded = s_B_excluded,
                  s_C_init = s_C_init, s_C_final = s_C_final)
  return(samples)
}

# units reallocated in each frame --> togliere e sostituire con riallocazione in B e in C
sample_scr_new_size_1 <- function(data, s_scr_reduced_size, n_A, n_B, n_C) {
  n <- n_A + n_B + n_C
  s_A_excluded <- s_scr_reduced_size$s_A_excluded
  s_A_final <- s_scr_reduced_size$s_A_final
  
  new_n_B <- n_B + nrow(s_A_excluded)/2
  s_B_init <- sample_n(data[data$B == 1 , ], new_n_B)
  s_B_final <- s_B_init %>% filter(domain == "b" | domain == "ab")
  s_B_excluded <- anti_join(s_B_init, s_B_final, by = c("id"))

  new_n_C <- n-(nrow(s_B_final) + nrow(s_A_final))
  s_C_init <- sample_n(data[data$C == 1 , ], new_n_C)
  s_C_final <- s_C_init 

  results <- list(s_A_excluded = s_A_excluded, s_B_excluded = s_B_excluded,
                   s_A_final = s_A_final, s_B_final = s_B_final, s_C_final = s_C_final)
  return(results)
}

# units reallocated only in less expensive frame, that here coincides with the last one
sample_scr_new_size_2 <- function(data, n_A, n_B, n_C) {
  n <- n_A + n_B + n_C
  s_A_init <- sample_n(data[data$A == 1 , ], n_A)
  s_A_final <- s_A_init %>% filter(domain == "a")
  s_A_excluded <- anti_join(s_A_init, s_A_final, by = c("id"))
  s_B_init <- sample_n(data[data$B == 1 , ], n_B)
  s_B_final <- s_B_init %>% filter(domain == "b" | domain == "ab")
  s_B_excluded <- anti_join(s_B_init, s_B_final, by = c("id"))
  n_C <- n - (nrow(s_A_final) + nrow(s_B_final))
  s_C_init <- sample_n(data[data$C == 1 , ], n_C)
  s_C_final <- s_C_init 
  samples <- list(s_A_init = s_A_init, s_A_final = s_A_final, s_A_excluded = s_A_excluded,
                  s_B_init = s_B_init, s_B_final = s_B_final, s_B_excluded = s_B_excluded,
                  s_C_init = s_C_init, s_C_final = s_C_final)
  return(samples)
}


# OVERLAP DESIGN
sample_mf <- function(data, n_A, n_B, n_C) {
  s_A <- sample_n(data[data$A == 1 , ], n_A)
  s_B <- sample_n(data[data$B == 1 , ], n_B)
  s_C <- sample_n(data[data$C == 1 , ], n_C)
  samples <- list(s_A_init = NULL, s_A_final = s_A, s_A_excluded = NULL,
                  s_B_init = NULL, s_B_final = s_B, s_B_excluded = NULL,
                  s_C_init = NULL, s_C_final = s_C)
  return(samples)
}


#### FUNZIONE DI COSTO

# contact_cost as a proportion of the full cost
cost <- function(s_scr_reduced_size, contact_cost = 0.1) {
  contact_A <- sum(s_scr_reduced_size$s_A_excluded$frame_cost_A*contact_cost)
  finished_A <- sum(s_scr_reduced_size$s_A_final$frame_cost_A)
  contact_B <- sum(s_scr_reduced_size$s_B_excluded$frame_cost_B*contact_cost)
  finished_B <- sum(s_scr_reduced_size$s_B_final$frame_cost_B)
  finished_C <- sum(s_scr_reduced_size$s_C_final$frame_cost_C)
  cost <- contact_A + finished_A + contact_B + finished_B + finished_C
  return(cost)
}


#### ESTIMATORS

# Screener design: stratified estimator
est_strat <- function(s, N_a, N_b_ab, N_C) {
  n_a <- nrow(s$s_A_final)
  n_b_ab <- nrow(s$s_B_final)
  n_C <- nrow(s$s_C_final) 
  hat_Y_a <- ifelse(n_a != 0, as.numeric((s$s_A_final %>% summarise(sum(Y)))), 0)
  hat_Y_b_ab <- ifelse(n_b_ab != 0, as.numeric((s$s_B_final %>% summarise(sum(Y)))), 0)
  hat_Y_C <- as.numeric((s$s_C_final %>% summarise(sum(Y))))
  hat_Y_str <- hat_Y_a*(N_a/n_a) + hat_Y_b_ab * (N_b_ab/n_b_ab) + hat_Y_C * (N_C/n_C)
  return(hat_Y_str)
}


# multiple frame design: multiplicity
est_mf_multiplicity <- function(s_mf, N_A, N_B, N_C) {
  hat_Y_m <- ((s_mf$s_A_final %>% filter( domain == "a") %>%  summarise(sum(Y))) +
                 (s_mf$s_A_final %>% filter( domain == "ab") %>% summarise(sum(Y)))/2 +
                 (s_mf$s_A_final %>% filter( domain == "ac") %>% summarise(sum(Y)))/2 +
                 (s_mf$s_A_final %>% filter( domain == "abc") %>% summarise(sum(Y)))/3)*(N_A/nrow(s_mf$s_A_final)) +
    ((s_mf$s_B_final %>% filter( domain == "b") %>% summarise(sum(Y))) +
       (s_mf$s_B_final %>% filter( domain == "ab") %>% summarise(sum(Y)))/2 +
       (s_mf$s_B_final %>% filter( domain == "bc") %>% summarise(sum(Y)))/2 +
       (s_mf$s_B_final %>% filter( domain == "abc") %>% summarise(sum(Y)))/3)*(N_B/nrow(s_mf$s_B_final)) +
    ((s_mf$s_C_final %>% filter( domain == "c") %>% summarise(sum(Y))) +
       (s_mf$s_C_final %>% filter( domain == "ac") %>% summarise(sum(Y)))/2 +
       (s_mf$s_C_final %>% filter( domain == "bc") %>% summarise(sum(Y)))/2 +
       (s_mf$s_C_final %>% filter( domain == "abc") %>% summarise(sum(Y)))/3)*(N_C/nrow(s_mf$s_C_final))
  
  return(hat_Y_m = as.numeric(hat_Y_m))
}

var_mf_multiplicity <- function(data, n_A, n_B, n_C) {
  data <- data %>% mutate(m = A + B + C)
  data_A <- data %>% filter(A == 1) 
  data_B <- data %>% filter(B == 1)
  data_C <- data %>% filter(C == 1)
  var_A <- (N_A - n_A)/(n_A*(N_A - 1))*(N_A*sum(data_A$Y^2/data_A$m^2) - (sum(data_A$Y/data_A$m))^2)
  var_B <- (N_B - n_B)/(n_B*(N_B - 1))*(N_B*sum(data_B$Y^2/data_B$m^2) - (sum(data_B$Y/data_B$m))^2)
  var_C <- (N_C - n_C)/(n_C*(N_C - 1))*(N_C*sum(data_C$Y^2/data_C$m^2) - (sum(data_C$Y/data_C$m))^2)
  var_m <- var_A + var_B + var_C
  return(var_m)
}


# multiple frame design: KA - coincide con multilicity quando probabilità di inclusione è n_q/N_q
est_mf_ka <- function(s_mf, N_A, N_B, N_C) {
  n_A <- nrow(s_mf$s_A_final)
  n_B <- nrow(s_mf$s_B_final)
  n_C <- nrow(s_mf$s_C_final)
  
  # single domains
  Y_a <- s_mf$s_A_final %>% filter(domain == "a") %>% summarise(sum(Y))
  Y_b <- s_mf$s_B_final %>% filter(domain == "b") %>% summarise(sum(Y))
  Y_c <- s_mf$s_C_final %>% filter(domain == "c") %>% summarise(sum(Y))
  # double domains
  Y_ab_A <- s_mf$s_A_final %>% filter(domain == "ab") %>% summarise(sum(Y))
  Y_ab_B <- s_mf$s_B_final %>% filter(domain == "ab") %>% summarise(sum(Y))
  Y_ac_A <- s_mf$s_A_final %>% filter(domain == "ac") %>% summarise(sum(Y))
  Y_ac_C <- s_mf$s_C_final %>% filter(domain == "ac") %>% summarise(sum(Y))
  Y_bc_B <- s_mf$s_B_final %>% filter(domain == "bc") %>% summarise(sum(Y))
  Y_bc_C <- s_mf$s_C_final %>% filter(domain == "bc") %>% summarise(sum(Y))
  # triple domains
  Y_abc_A <- s_mf$s_A_final %>% filter(domain == "abc") %>% summarise(sum(Y))
  Y_abc_B <- s_mf$s_B_final %>% filter(domain == "abc") %>% summarise(sum(Y))
  Y_abc_C <- s_mf$s_C_final %>% filter(domain == "abc") %>% summarise(sum(Y))
  # KA
  hat_Y_ka <- 
    (Y_a + Y_ab_A*((n_A/N_A)/(n_A/N_A + n_B/N_B)) + Y_ac_A*((n_A/N_A)/(n_A/N_A + n_C/N_C)) + 
       Y_abc_A*((n_A/N_A)/(n_A/N_A + n_B/N_B + n_C/N_C)))*(N_A/n_A) +
    (Y_b + Y_ab_B*((n_B/N_B)/(n_A/N_A + n_B/N_B)) + Y_bc_B*((n_B/N_B)/(n_B/N_B + n_C/N_C)) 
     + Y_abc_B*((n_B/N_B)/(n_A/N_A + n_B/N_B + n_C/N_C)))*(N_B/n_B) +  
    (Y_c + Y_bc_C*((n_C/N_C)/(n_B/N_B + n_C/N_C)) + Y_ac_C*((n_C/N_C)/(n_A/N_A + n_C/N_C)) + 
       Y_abc_C*((n_C/N_C)/(n_A/N_A + n_B/N_B + n_C/N_C)))*(N_C/n_C)
  
  return(hat_Y_ka = as.numeric(hat_Y_ka))
}


var_mf_ka <- function(data, n_A, n_B, n_C) {
}



### CALCULATING SAMPLE SIZE


## proportional allocation in stratified/screener design

proportional_n_screener <- function(N, N_a, N_b_ab, N_C, n) {
  n_a <- ceiling((N_a/N)*n)
  n_b_ab <- ceiling((N_b_ab/N)*n)
  n_C <- n - (n_a + n_b_ab)
  return(list(n_a = n_a, n_b_ab = n_b_ab, n_C = n_C))
}


## optimal allocation in stratified/screener design
# first we need to calculate variances
sigma2_strata <- function(data, N_A, N_B, N_C) {
  data_a <- data %>% filter(A == 1 & B == 0 & C == 0) 
  data_ab_b <- data %>% filter((A == 1 & B == 1 & C == 0) | (A == 0 & B == 1 & C == 0))
  data_C <- data %>% filter(C == 1)
  sigma2_a <- var(data_a$Y)
  sigma2_ab_b <- var(data_ab_b$Y)
  sigma2_C <- var(data_C$Y)
  return(sigma2 = list(sigma2_a = sigma2_a, sigma2_ab_b = sigma2_ab_b,
                       sigma2_C = sigma2_C))
}

# second we can calculate the optimal sample size
optimal_n_screener <- function(N, N_a, N_A, N_b_ab, N_B, N_C, c_a, c_b_ab, c_C, C, sigma2) {
  V_a <- N_a*N_A*sigma2[["sigma2_a"]]
  V_b_ab <- N_b_ab*N_B*sigma2[["sigma2_ab_b"]]
  V_C <- N_C^2*sigma2[["sigma2_C"]]
  n_a <- C*(sqrt(V_a/c_a)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C))) 
  n_b_ab <- C*(sqrt(V_b_ab/c_b_ab)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C)))
  n_C <- C*(sqrt(V_C/c_C)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C)))
  return(list(n_a = n_a, n_b_ab = n_b_ab, n_C = n_C))
}


## proportional allocation in mf

proportional_n_mf <- function(N, N_A, N_B, N_C, n) {
  n_A <- ceiling((N_A/(N_A+N_B+N_C))*n)
  n_B <- ceiling((N_B/(N_A+N_B+N_C))*n)
  n_C <- n - (n_A + n_B)
  return(list(n_A = n_A, n_B = n_B, n_C = n_C))
}


## optimal sample size allocation (with respect to variance)

# first, we need multiplicity adjusted frame variance: I am calculating exact vars
sigma_alpha_values <- function(data, N_A, N_B, N_C) {
  data <- data %>% mutate(m = A + B + C)
  data_A <- data %>% filter(A == 1) 
  data_B <- data %>% filter(B == 1)
  data_C <- data %>% filter(C == 1)
  sigma2_alpha_A <- mean(data_A$Y^2/data_A$m^2) - (mean(data_A$Y/data_A$m))^2
  sigma2_alpha_B <- mean(data_B$Y^2/data_B$m^2) - (mean(data_B$Y/data_B$m))^2
  sigma2_alpha_C <- mean(data_C$Y^2/data_C$m^2) - (mean(data_C$Y/data_C$m))^2
  return(sigma2_alpha = list(sigma2_alpha_A = sigma2_alpha_A, sigma2_alpha_B = sigma2_alpha_B,
                            sigma2_alpha_C = sigma2_alpha_C))
}


# second, we calculate the optimal sample size for each frame: equation (22) corretta, draft Mecatti
optimal_n_mf <- function(N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha) {
  
  n_q_opt_denominator <- (N_A*sqrt((N_A/(N_A-1))*sigma2_alpha[["sigma2_alpha_A"]]*cost_frames[1])) +
                         (N_B*sqrt((N_B/(N_B-1))*sigma2_alpha[["sigma2_alpha_B"]]*cost_frames[2])) +
                         (N_C*sqrt((N_C/(N_C-1))*sigma2_alpha[["sigma2_alpha_C"]]*cost_frames[3]))
  
  # (C_c_0) that multiplies everything actually translates into the total sample size n... è necessario? non è meglio
  # metterci n totale allora?
  n_A_opt <- ((C-c_0)*(N_A*sqrt((N_A/(N_A-1))*sigma2_alpha[["sigma2_alpha_A"]]))/sqrt(cost_frames[1]))/n_q_opt_denominator
  n_B_opt <- ((C-c_0)*(N_B*sqrt((N_B/(N_B-1))*sigma2_alpha[["sigma2_alpha_B"]]))/sqrt(cost_frames[2]))/n_q_opt_denominator
  n_C_opt <- ((C-c_0)*(N_C*sqrt((N_C/(N_C-1))*sigma2_alpha[["sigma2_alpha_C"]]))/sqrt(cost_frames[3]))/n_q_opt_denominator
  n <- n_A_opt + n_B_opt + n_C_opt
  
  n_results <- list(n_A = n_A_opt, n_B = n_B_opt, n_C = n_C_opt, n = n)
  return(n_results)
}  


n_size <- function(type = "proportional", N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha, n) {
  switch(type,
         proportional = proportional_n_mf(N, N_A, N_B, N_C, n),
         optimal_eq_mecatti = optimal_n_mf(N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha),
         # non sono sicura sia utile implementarla, al momento non l'ho implementata
         optimal_eq_lohr = optimal_n_mf2(N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha)) 
}

