require(dplyr)

### POPULATION

# generate pop with different variance differntials
# at the moment, N fixed and overlap fixed!!

generate_pop <- function(N, p_domains, mu_domains, sd_domains, frame_cost) {
  
  p_domains = p_domains/sum(p_domains)
  
  n_domains = round(N*p_domains)
  
  domains = rbind(matrix(rep(c(1,0,0), n_domains[1]), ncol=3, byrow=TRUE),
                  matrix(rep(c(0,1,0), n_domains[2]), ncol=3, byrow=TRUE),
                  matrix(rep(c(0,0,1), n_domains[3]), ncol=3, byrow=TRUE),
                  matrix(rep(c(1,1,0), n_domains[4]), ncol=3, byrow=TRUE),
                  matrix(rep(c(0,1,1), n_domains[5]), ncol=3, byrow=TRUE),
                  matrix(rep(c(1,0,1), n_domains[6]), ncol=3, byrow=TRUE),
                  matrix(rep(c(1,1,1), n_domains[7]), ncol=3, byrow=TRUE))
  
  colnames(domains) = c("A","B","C")
  
  population <- data.frame(
    id = 1:N,
    domains
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
  
  N_a   = n_domains[1]
  N_b   = n_domains[2]
  N_c   = n_domains[3]
  N_ab  = n_domains[4]
  N_bc  = n_domains[5]
  N_ac  = n_domains[6]
  N_abc = n_domains[7]
  
  set.seed(123)
  
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


#N <- 70

#mu_domains = rep(0,7)
#sd_domains = rep(1,7)
#frame_cost = 1:3

#p_domains = c(2,2,2,1,1,1,1) # a, b, c, ab, bc, ac, abc

#generate_pop(N, p_domains, mu_domains, sd_domains, frame_cost)


### FINE TUNE DOMAIN VARIANCE

# Is the formula for the variance of 
# a single frame based on the domains

formula_var <- function(pop_d, mu_domains, sd_domains){
    pop_d <- pop_d/sum(pop_d)
  result <- pop_d%*%(sd_domains^2) + 
    (pop_d%*%(mu_domains^2) - (pop_d%*%mu_domains)^2)
  return(as.numeric(result))
}

# Domain order:
# a - b - c - ab - bc - ac - abc

frame_sd <- function(pop_d, mu_domains, sd_domains){
  pos_select <- list(c(1, 4, 6:7), # Frame A
                     c(2, 4:5, 7), # Frame B
                     c(3, 5:7)) # Frame C
  lapply(pos_select, 
        function(x) sqrt(formula_var(pop_d = pop_d[x],
                                 mu_domains = mu_domains[x],
                                 sd_domains = sd_domains[x])))
}

# Toy example to fine tune the sd of each frame
# n_doms <- c(rep(1e3, 3), 1e2, 2e3, 4e2, 1e2) # Population of each domain
# mu_s <- c(rep(5, 3), 4, 3, 1, 1) # Means
# sd_s <- c(rep(9, 3), 4, 5, 6, 6) # Standard deviations
# frame_sd(pop_d = n_doms, mu_domains = mu_s, sd_domains = sd_s)
  
### SCREENER DESIGN: implementation of sampling design

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
  s_C_excluded <- anti_join(s_C_init, s_C_final, by = c("id"))

  samples <- list(s_A_init = s_A_init, s_A_final = s_A_final, s_A_excluded = s_A_excluded,
                  s_B_init = s_B_init, s_B_final = s_B_final, s_B_excluded = s_B_excluded,
                  s_C_init = s_C_init, s_C_final = s_C_final, s_C_excluded = s_C_excluded)
  
  return(samples)
}

#s_scr_reduced_size <- sample_scr_reduced_size(data = population, n_A = n_A, n_B = n_B, n_C = n_C)

# method 1: allocate in B and C
sample_scr_new_size_1 <- function(data, n_A, n_B, n_C) {
  n <- n_A + n_B + n_C

  s_A_init <- sample_n(data[data$A == 1 , ], n_A)
  s_A_final <- s_A_init %>% filter(domain == "a")
  s_A_excluded <- anti_join(s_A_init, s_A_final, by = c("id"))
  
  n_B_plus <- n_B + (nrow(s_A_excluded))
  s_B_init <- sample_n(data[data$B == 1 , ], n_B_plus)
  s_B_final <- s_B_init %>% filter(domain == "b" | domain == "ab")
  s_B_excluded <- anti_join(s_B_init, s_B_final, by = c("id"))
  
  n_C_plus <- n_C + (nrow(s_B_excluded))
  s_C_init <- sample_n(data[data$C == 1 , ], n_C_plus)
  s_C_final <- s_C_init
  s_C_excluded <- anti_join(s_C_init, s_C_final, by = c("id"))
  
  samples <- list(s_A_init = s_A_init, s_A_final = s_A_final, s_A_excluded = s_A_excluded,
                  s_B_init = s_B_init, s_B_final = s_B_final, s_B_excluded = s_B_excluded,
                  s_C_init = s_C_init, s_C_final = s_C_final, s_C_excluded = s_C_excluded)
  
  return(samples)
}

# method 2: units reallocated only in less expensive frame, that here coincides with the last one
sample_scr_new_size_2 <- function(data, n_A, n_B, n_C) {
  n <- n_A + n_B + n_C
  
  s_A_init <- sample_n(data[data$A == 1 , ], n_A)
  s_A_final <- s_A_init %>% filter(domain == "a")
  s_A_excluded <- anti_join(s_A_init, s_A_final, by = c("id"))
  
  s_B_init <- sample_n(data[data$B == 1 , ], n_B)
  s_B_final <- s_B_init %>% filter(domain == "b" | domain == "ab")
  s_B_excluded <- anti_join(s_B_init, s_B_final, by = c("id"))
  
  n_C_plus <- n - (nrow(s_A_final) + nrow(s_B_final))
  s_C_init <- sample_n(data[data$C == 1 , ], n_C_plus)
  s_C_final <- s_C_init
  s_C_excluded <- anti_join(s_C_init, s_C_final, by = c("id"))
  
  samples <- list(s_A_init = s_A_init, s_A_final = s_A_final, s_A_excluded = s_A_excluded,
                  s_B_init = s_B_init, s_B_final = s_B_final, s_B_excluded = s_B_excluded,
                  s_C_init = s_C_init, s_C_final = s_C_final, s_C_excluded = s_C_excluded)
  return(samples)
}


# OVERLAP DESIGN
sample_mf <- function(data, n_A, n_B, n_C) {
  s_A <- sample_n(data[data$A == 1 , ], n_A)
  s_B <- sample_n(data[data$B == 1 , ], n_B)
  s_C <- sample_n(data[data$C == 1 , ], n_C)
  samples <- list(s_A_init = NULL, s_A_final = s_A, s_A_excluded = NULL,
                  s_B_init = NULL, s_B_final = s_B, s_B_excluded = NULL,
                  s_C_init = NULL, s_C_final = s_C, s_C_excluded = NULL)
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


var_strat <- function(data, s) {
  data_a <- data %>% filter(domain == "a")
  sigma2_a <- ifelse(nrow(data_a) != 0, var(data_a$Y)*((nrow(data_a)-1)/nrow(data_a)), 0)
  N_a <- nrow(data_a)
  n_a <- nrow(s$s_A_final)
  
  data_b_ab <- data %>% filter(domain == "a" | domain == "ab")
  sigma2_b_ab <- ifelse(nrow(data_b_ab) != 0 , var(data_b_ab$Y)*((nrow(data_b_ab)-1)/nrow(data_b_ab)), 0)
  N_b_ab <- nrow(data_b_ab)
  n_b_ab <- nrow(s$s_B_final)
  
  data_C <- data %>% filter(C == 1)
  sigma2_C <- ifelse(nrow(data_C) != 0, var(data_C$Y)*((nrow(data_C)-1)/nrow(data_C)), 0)
  N_C <- nrow(data_C)
  n_C <- nrow(s$s_C_final)
  
  var_Y_str <- (N_a^2*(sigma2_a/n_a) + N_b_ab^2*(sigma2_b_ab/n_b_ab) + N_C^2*(sigma2_C/n_C)) -
    (N_a*sigma2_a + N_b_ab*sigma2_b_ab + N_C*sigma2_C)
  
  return(var_Y_str)
}


# multiple frame design: multiplicity estimator
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

# exact variance of multiplicity estimator
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


# multiple frame design: KA estimator - coincide con multilicity quando probabilità di inclusione è n_q/N_q
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


# exact variance of kalton anderson estimator: NOT YET IMPLEMENTED
var_mf_ka <- function(data, n_A, n_B, n_C) {
}


# Hartley estimator for three frames: NOT YET IMPLEMENTED
# to be implemented following equation 18 and 17 in Mecatti and Singh 2014
# est_mf_h


### SAMPLE SIZE CALCULATION


## 1) PROPORTIONAL
## proportional allocation in stratified/screener design
proportional_n_screener <- function(N, N_a, N_b_ab, N_C, n) {
  n_a <- ceiling((N_a/N)*n)
  n_b_ab <- ceiling((N_b_ab/N)*n)
  n_C <- n - (n_a + n_b_ab)
  return(list(n_a = n_a, n_b_ab = n_b_ab, n_C = n_C))
}


## proportional allocation in mf
proportional_n_mf <- function(N, N_A, N_B, N_C, n) {
  n_A <- ceiling((N_A/(N_A+N_B+N_C))*n)
  n_B <- ceiling((N_B/(N_A+N_B+N_C))*n)
  n_C <- n - (n_A + n_B)
  return(list(n_A = n_A, n_B = n_B, n_C = n_C))
}


## OPTIMAL SAMPLE SIZE ALLOCATION (MINIMIZE ESTIMATOR VARIABILITY)

## optimal allocation in stratified/screener design
# first we need to calculate variances: this function is only useful for finding
# strata variances to be included in the optimal allocation calculation
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
# this function returns the sample size to be allocated in each strata
optimal_n_screener <- function(N, N_a, N_A, N_b_ab, N_B, N_C, c_a, c_b_ab, c_C, C, sigma2) {
  V_a <- N_a*N_A*sigma2[["sigma2_a"]]
  V_b_ab <- N_b_ab*N_B*sigma2[["sigma2_ab_b"]]
  V_C <- N_C^2*sigma2[["sigma2_C"]]
  n_a <- C*(sqrt(V_a/c_a)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C))) 
  n_b_ab <- C*(sqrt(V_b_ab/c_b_ab)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C)))
  n_C <- C*(sqrt(V_C/c_C)/(sqrt(V_a*c_a) + sqrt(V_b_ab*c_b_ab) + sqrt(V_C/c_C)))
  return(list(n_a = n_a, n_b_ab = n_b_ab, n_C = n_C))
}

## optimal sample size allocation in overlap design
# first, we need multiplicity adjusted frame variance: I am calculating exact vars
# this function is only useful for finding frames variances to be included in the optimal allocation calculation

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


# if we choose type = proportional, the functions that calculate proprtional allocation are used
# if we choose optimal allocation, the functions for optimal allocation are used
n_size <- function(type = "proportional", N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha, n) {
  switch(
    type,
    proportional = proportional_n_mf(N, N_A, N_B, N_C, n),
    optimal = optimal_n_mf(N, N_A, N_B, N_C, C, c_0, cost_frames, sigma2_alpha)
    )
}

