require(dplyr)

### POPULATION

# generate pop with different variance differntials
# at the moment, N fixed and overlap fixed!!

generate_pop <- function(shape_pars_frame, frame_cost) {
  N <- 70000
  
  A <- c(rep(1, 40000), rep(0, 30000))
  B <- c(rep(0, 10000), rep(1, 10000), rep(0, 10000), rep(1, 30000), rep(0, 10000))
  C <- c(rep(0, 20000), rep(1, 30000), rep(0, 10000), rep(1, 10000))
  
  population <- data.frame(
    id = 1:N,
    A = A,
    B = B,
    C = C,
    Y = ifelse( A == 1, rgamma(10, shape_pars_frame[1]), 
                               ifelse(B == 1, rgamma(10, shape_pars_frame[2]),
                                                      ifelse(C == 1, rgamma(10, shape_pars_frame[3]), 0))) 
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


# units reallocated in each frame
sample_scr_new_size_1 <- function(data, s_scr_reduced_size, n_A, n_B, n_C) {
  s_A_excluded <- s_scr_reduced_size$s_A_excluded
  s_A_final <- s_scr_reduced_size$s_A_final
  while(nrow(s_A_final) < n_A) {
    new_n_A <- n_A - nrow(s_A_final)
    s_A_init_new <- sample_n(data[data$A == 1 , ] %>% anti_join(s_A_final), new_n_A)
    s_A_final_new <- s_A_init_new %>% filter(domain == "a")
    s_A_excluded <- s_A_excluded %>% bind_rows(anti_join(s_A_init_new, s_A_final_new))
    s_A_final <- s_A_final %>% bind_rows(s_A_final_new)
  }
  
  s_B_excluded <- s_scr_reduced_size$s_B_excluded
  s_B_final <- s_scr_reduced_size$s_B_final
  while(nrow(s_B_final) < n_B) {
    new_n_B <- n_B - nrow(s_B_final)
    s_B_init_new <- sample_n(data[data$B == 1 , ] %>% anti_join(s_B_final), new_n_B)
    s_B_final_new <- s_B_init_new %>% filter(domain == "b" | domain == "ab")
    s_B_excluded <- s_B_excluded %>% bind_rows(anti_join(s_B_init_new, s_B_final_new))
    s_B_final <- s_B_final %>% bind_rows(s_B_final_new)
  }
  
  s_C_init <- s_scr_reduced_size$s_C_init
  s_C_final <- s_scr_reduced_size$s_C_final
  
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


# MULTIPLE FRAME DESIGN
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
est_strat <- function(s) {
  n_A <- nrow(s$s_A_final)
  n_B <- nrow(s$s_B_final)
  n_C <- nrow(s$s_C_final) 
  hat_Y_A <- ifelse(n_A != 0, as.numeric((s$s_A_final %>% summarise(sum(Y)))/n_A), 0)
  hat_Y_B <- ifelse(n_B != 0, as.numeric((s$s_B_final %>% summarise(sum(Y)))/n_B), 0)
  hat_Y_C <- as.numeric((s$s_C_final %>% summarise(sum(Y)))/n_C)
  hat_Y_str <- hat_Y_A*1000/7000 + hat_Y_B * 2000/7000 + hat_Y_C * 4000/7000
  return(hat_Y_str)
}


# multiple frame design: multiplicity
est_mf_multiplicity <- function(s_mf) {
  hat_Y_m <- ((s_mf$s_A_final %>% filter( domain == "a") %>%  summarise(sum(Y))) +
                 (s_mf$s_A_final %>% filter( domain == "ab") %>% summarise(sum(Y)))/2 +
                 (s_mf$s_A_final %>% filter( domain == "ac") %>% summarise(sum(Y)))/2 +
                 (s_mf$s_A_final %>% filter( domain == "abc") %>% summarise(sum(Y)))/3)/nrow(s_mf$s_A_final) +
    ((s_mf$s_B_final %>% filter( domain == "b") %>% summarise(sum(Y))) +
       (s_mf$s_B_final %>% filter( domain == "ab") %>% summarise(sum(Y)))/2 +
       (s_mf$s_B_final %>% filter( domain == "bc") %>% summarise(sum(Y)))/2 +
       (s_mf$s_B_final %>% filter( domain == "abc") %>% summarise(sum(Y))/3))/nrow(s_mf$s_B_final) +
    ((s_mf$s_C_final %>% filter( domain == "c") %>% summarise(sum(Y))) +
       (s_mf$s_C_final %>% filter( domain == "ac") %>% summarise(sum(Y)))/2 +
       (s_mf$s_C_final %>% filter( domain == "bc") %>% summarise(sum(Y)))/2 +
       (s_mf$s_C_final %>% filter( domain == "abc") %>% summarise(sum(Y)))/3)/nrow(s_mf$s_C_final)
  
  return(hat_Y_m = as.numeric(hat_Y_m))
}


# multiple frame design: KA - coincide con multilicity per pi n_A/N_A
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
    ((Y_a + Y_ab_A*((n_A/N_A)/(n_A/N_A + n_B/N_B)) + Y_ac_A*((n_A/N_A)/(n_A/N_A + n_C/N_C)) + 
       Y_abc_A*((n_A/N_A)/(n_A/N_A + n_B/N_B + n_C/N_C)))/n_A +
    (Y_b + Y_ab_B*((n_B/N_B)/(n_A/N_A + n_B/N_B)) + Y_bc_B*((n_B/N_B)/(n_B/N_B + n_C/N_C)) 
     + Y_abc_B*((n_B/N_B)/(n_A/N_A + n_B/N_B + n_C/N_C)))/n_B +  
    (Y_c + Y_bc_C*((n_C/N_C)/(n_B/N_B + n_C/N_C)) + Y_ac_C*((n_C/N_C)/(n_A/N_A + n_C/N_C)) + 
       Y_abc_C*((n_C/N_C)/(n_A/N_A + n_B/N_B + n_C/N_C)))/n_C)
  
  return(hat_Y_ka = as.numeric(hat_Y_ka))
}



### CALCULATING SAMPLE SIZE


## proportional : stratified

proportional_n_screener <- function(N, N_a, N_b_ab, N_C, n) {
  n_a <- ceiling((N_a/N)*n)
  n_b_ab <- ceiling((N_b_ab/N)*n)
  n_C <- n - (n_a + n_b_ab)
  return(list(n_a = n_a, n_b_ab = n_b_ab, n_C = n_C))
}

## proportional : MF

proportional_n_mf <- function(N, N_A, N_B, N_C, n) {
  n_A <- ceiling((N_A/(N_A+N_B+N_C))*n)
  n_B <- ceiling((N_B/(N_A+N_B+N_C))*n)
  n_C <- n - (n_A + n_B)
  return(list(n_A = n_A, n_B = n_B, n_C = n_C))
}


# optimal (with respect to variance)

