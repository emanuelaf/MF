require(dplyr)

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

# how do I find lambda?
optimal_sample_size <- function(N_A, N_a, lambda, N_ab, ov_c_A, N_B, N_b, ov_c_B) {
  alpha <- sqrt((N_A*(N_a+lambda^2*N_ab)/ov_c_A))*(sqrt(N_A*(N_a + lambda^2*N_ab)/ov_c_A) + 
                                                     sqrt(N_B*(N_b+lambda^2*N_ab)/ov_c_B))
  return(list(alpha = alpha))
}

optimal_sample_size(N_A = 10, N_a = 1, lambda = 2, N_ab = 3, ov_c_A = 3, N_B = 4, N_b = 1, ov_c_B = 2)

#######################
### overlap design ####
#######################

# how does the overlap happen?
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

ov_variance_A <- function(N_a, N_A, S_a, lambda, N_ab, S_ab) {
  variance_A <- N_a*N_A*S_a^2 + lambda^2*N_ab*N_A*S_ab^2
  return(list(variance_A = variance_A))
}

ov_variance_B <- function(N_n, N_B, S_b, lambda, N_ab, S_ab) {
  variance_B <- N_b*N_B*S_b^2 + (1-lambda)^2*N_ab*N_B*S_ab^2
  return(list(variance_B = variance_B))
}

# how to calculate hartley's variance?
variance_hartley <- function() {
  
}


#########################
### screener design #####
#########################

scr_sampling <- function(data, n, n_B) {
  n_A <- n - n_B
  s_B <- sample_n(data[data$B == 1 , ], n_B, replace = FALSE)
  s_b <- s_B %>% filter(domain == "b")
  s_A <- sample_n(data[data$A == 1 , ], n_A, replace = FALSE)
  final_sample <- list(s_A = s_A, s_b = s_b)
  return(final_sample)
}

scr_sampling(data = data, n = 10, n_B = 5)

