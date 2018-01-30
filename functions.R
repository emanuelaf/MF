require(dplyr)

# given quantities:
# N, N_A, N_B, N_a, N_ab, N_b, n
# f_a = N_a/N, f_b = N_b/N

##############################
##### set up functions #######
##############################

produce_dataframe <- function(N, f_a, f_b, overlap_size, n_frames) {
  y <- rnorm(N)
  id <- 1:N
  if (n_frames == 2) {
    N_A <- N*f_a+N*overlap_size
    N_B <- N*f_b+N*overlap_size
    id_A <- sort(sample(1:N, N_A))
    id_B <- sort(c(id[!(id %in% id_A)], sample(id_A, N*overlap_size)))
    df <- data.frame(id = 1:N, y = y, A = numeric(N), B = numeric(N))
    df[id_A,]$A <- 1 
    df[id_B,]$B <- 1
    df$domain[df$A == 1 & df$B == 0] <-  "a" 
    df$domain[df$A == 0 & df$B == 1] <-  "b" 
    df$domain[df$A == 1 & df$B == 1] <-  "ab"
    return(df)
  }
  if (n_frames != 2) {
    print("non ancora sviluppato quando numero di frames > 2")
  }

}

data <- produce_dataframe(N = 100, f_a = 0.09, f_b = 0.4, overlap_size = 0.51, n_frame = 2)

# still to do: function to optimise lambda



# 

# sample size to be allocated in A and B given n
optimal_sample_size <- function(N_A, N_a, n, lambda = "lambda_opt", N_ab, ov_c_A, N_B, N_b, ov_c_B) {
  if (lambda == "lambda_opt") {print("to add function for lambda optimisation")}
  alpha <- (sqrt(N_A*(N_a+lambda^2*N_ab)/ov_c_A))/((sqrt(N_A*(N_a+lambda^2*N_ab)/ov_c_A)) + 
                                                     (sqrt(N_B*(N_b+(1-lambda)^2*N_ab)/ov_c_B)))
  return(list(n_A = alpha*n, n_B = (1-alpha)*n, alpha = alpha))
}

optimal_sample_size(N_A = 10, N_a = 5, lambda = 0.5, n = 4, 
                    N_ab = 10, ov_c_A = 3, N_B = 10, N_b = 5, ov_c_B = 2)

#######################
####### designs #######
#######################

# overlap, how does the overlap happen?
ov_sampling <- function(data, n, n_A) {
  n_B <- n - n_A
  s_A <- sample_n(data[data$A == 1 , ], n_A, replace = FALSE)
  s_B <- sample_n(data[data$B == 1 , ], n_B, replace = FALSE)
  s_a <- s_A %>% filter(domain == "a")
  s_ab <- s_A %>% filter(domain == "ab") %>% bind_rows(s_B %>% filter(domain == "ab"))
  s_b <- s_B %>% filter(domain == "b")
  final_sample <- list(s_a = s_a, s_b = s_b, s_ab = s_ab)
  return(final_sample)
}

ov_sampling(data = data, n = 10, n_A = 5)


# screening
scr_sampling <- function(data, n, n_B) {
  n_A <- n - n_B
  s_B <- sample_n(data[data$B == 1 , ], n_B, replace = FALSE)
  s_b <- s_B %>% filter(domain == "b")
  s_A <- sample_n(data[data$A == 1 , ], n_A, replace = FALSE)
  final_sample <- list(s_A = s_A, s_b = s_b)
  return(final_sample)
}

scr_sampling(data = data, n = 10, n_B = 5)


###########################
### estimation ####
##########################

#### overlap
# n_a/N_a because we consider a simple random sampling (following Lohr and Brick)
est_hat_Y_a <-  function(s, n_A, N_A) {
  s_a <- s %>% filter(domain == "a")
  hat_Y_a <- sum(s_a$Y)*n_a/N_a
  return(hat_Y_a)
}

# this is no good because I need to recognise who came from sample A and who from sample B
est_hat_Y_ab_A <- function(s, N_ab_A) {
  s_ab <- s %>% filter(domain == "ab")
  hat_Y_ab_A <- sum(s_a$Y)*n_a/N_a
  return(hat_Y_ab_A)
} 

# this is no good because I need to recognise who came from sample A and who from sample B
est_hat_Y_ab_B <- function(s, N_ab_B) {
  s_ab <- s %>% filter(domain == "ab")
  hat_Y_ab_B <- sum(s_a$Y)*n_a/N_a
  return(hat_Y_ab_B)
} 

est_hat_Y_b <- function(s, n_B, N_B) {
  s_b <- s %>% filter(domain == "b")
  hat_Y_b <- sum(s_b$Y)*n_b/N_b
  return(hat_Y_b)
} 


# overall estimation
hat_Y_lambda <- function(hat_Y_a, hat_Y_ab_A, hat_Y_ab_B, hat_Y_b, lambda) {
  hat_Y_lambda_est <- hat_Y_a + lambda*hat_Y_ab_A + (1-lambda)*hat_Y_ab_B + hat_Y_b
}

# variance estimation S_a, S_ab and S_b, where do they come from??
ov_variance_hat_Y_A <- function(N_a, N_A, S_a, lambda, N_ab, S_ab) {
  var_hat_Y_A <- N_a*N_A*S_a^2 + lambda^2*N_ab*N_A*S_ab^2
  return(list(var_hat_Y_A = var_hat_Y_A))
}

ov_variance_hat_Y_B <- function(N_n, N_B, S_b, lambda, N_ab, S_ab) {
  var_hat_Y_B <- N_b*N_B*S_b^2 + (1-lambda)^2*N_ab*N_B*S_ab^2
  return(list(var_hat_Y_B = var_hat_Y_B))
}


variance_hat_Y_lambda <- function(n_A, n_B, variance_A, variance_B) {
  var_hat_Y_lambda <- variance_A/n_A + variance_B/n_B
  return(list(var_hat_Y_lambda = var_hat_Y_lambda))
}



### screening
# n_a/N_a because we consider a simple random sampling (following Lohr and Brick)
scr_est_hat_Y_A <-  function(s, n_A, N_A) {
  s_A <- s %>% filter_if(str_detect(domain, 'a'))
  scr_hat_Y_A <- sum(s_A$Y)*n_A/N_A
  return(scr_hat_Y_A)
}


scr_est_hat_Y_b <-  function(s, n_b, N_b) {
  s_b <- s %>% filter_if(domain == 'b')
  scr_hat_Y_b <- sum(s_b$Y)*n_b/N_b
  return(scr_hat_Y_b)
}

# overall estimation
scr_est_hat_Y <-  function(scr_hatY_A, scr_hat_Y_b) {
  scr_hat_Y <- scr_hatY_A + scr_hat_Y_b
  return(scr_hat_Y)
}


scr_variance_hat_Y_A  <- function(N_A, S_A, n_A) {
  scr_var_hat_Y_A <- (N_A^2*S_A^2)/n_A
  return(scr_var_hat_Y_A)
}


scr_variance_hat_Y_b  <- function(N_b, S_b, n_b) {
  scr_var_hat_Y_b <- (N_b^2*S_b^2)/n_b
  return(scr_var_hat_Y_b)
}
