##################################
#
# Multivariate time warping
#
# 30/09/19
#
##################################

library(tidyverse)
library(minqa)
source("parscale.R")

# Important Functions -----------------------------------------------------

# Function defining basis splines

B = function(x, j, n, K) {
  
  # Creating j-th spline
  
  b <- numeric(length(x))
  
  b[x < K[j]] <- 0
  b[x >= K[j] & x <= K[j + 1]] <- (x[x >= K[j] & x <= K[j + 1]] - K[j]) / (K[j + 1] - K[j])
  b[x > K[j + 1]] <- 1
  
  return(b)
  
}

# Defining warping function phi

phi <- function(s, a, order = 1) {
  
  order = 1
  
  # Creating knots vector 
  
  n_knots <- length(a) - 1
  knots <- seq(1, by = round(length(s) / n_knots), length.out = n_knots)
  knots <- c(knots, length(s))
  
  # Padding knots vector to make right length
  
  req_length <- n_knots + order + 1
  zero <- FALSE
  while (length(knots) < req_length) {
    
    if (zero) {
      knots <- c(1, knots)
    } else {
      knots <- c(knots, length(s))
    }
    
    zero <- !zero
    
  }
  
  # Creating t_hat
  
  t_hat <- a[1] * B(s, 1, order, c(0, knots))
  
  for (i in 1:n_knots) {
    
    t_hat <- t_hat + a[i + 1] * B(s, i, order, knots)
    
  }
  
  return(list(
    t_hat = t_hat,
    knots = knots
  ))
  
}

# Function to get closest depth from t and assign it to t_hat

match_depths <- function(t_hat, t) {
  
  for (i in 1:length(t_hat)) {
    
    distances <- abs(t_hat[i] - t)
    index <- which.min(distances)
    t_hat[i] <- t[index]
    
  } 
  
  return(t_hat)
  
}

# Function for finding correlation between X(s) and Y(phi(s)), with a penalty for optimisation

corr_mv <- function(theta, s, t, lambda_1, lambda_2, lambda_3) {
  
  a <- theta
  
  t_hat <- phi(s, a, order)$t_hat
  t_hat_matched <- match_depths(t_hat, t)
  transformation <- sapply(t_hat_matched, function(x) {which(t == x)})
  
  difference_penalty <- sum(diff(a[-1])^2)
  decrease_penalty <- sum(transformation %>% diff() <= 0)
  
  # Penalty for major deviation from average mm per year
  slope_penalty <- sum((a / 12 - 0.677919 / 1000)^2)
  
  Y_hat <- data[transformation, ]
  # Multiply by -1 to create minimisation problem
  (sum(abs(cor(sunspots, Y_hat))) - lambda_1 * difference_penalty - lambda_2 * slope_penalty - lambda_3 * decrease_penalty) * -1
  
}

# Function for doing the time warp

time_warp_mv <- function(a_0, s, t, sunspots, data, parscale, lambda_1 = 0, lambda_2 = 0, lambda_3 = 0, i = 1, plots = TRUE) {
  
  order <- 1
  
  # Optimisation of a using optim
  result <- optim(a_0, corr_mv, s = s, t = t, lambda_1 = lambda_1, lambda_2 = lambda_2, lambda_3 = lambda_3, control = list(maxit = 1000, parscale = parscale))
  
  # Results
  
  a <- result$par
  
  spline <- phi(s, a, order) 
  t_hat <- spline$t_hat
  t_hat_matched <- match_depths(t_hat, t)
  transformation <- sapply(t_hat_matched, function(x) {which(t == x)})
  Y_hat <- data[transformation, ]
  t_warped <- t[transformation]
  
  if (plots) {
    
    # Plots
    
    par(mfrow = c(3, 2))
    
    # Plot 1
    plot(sunspots, Y_hat[, i], ylab = bquote(hat(Y)(t)), xlab = "X(s)")
    
    # Plot 2
    plot(t, data[, i], type = "l", ylab = "Y(t)")
    
    # Plot 3
    plot(s, t_hat, type = "l", ylab = bquote(hat(phi)(s)))
    abline(v = spline$knots, lty = 2)
    
    # Plot 4
    plot(s, sunspots, type = "l", ylab = "X(s)")
    
    # Plot 5
    plot(1:length(transformation), t[transformation], type = "l", xlab = "s", ylab = bquote(hat(t)))
    lines(1:length(transformation), min_depth + (max_depth - min_depth) / length(sunspots) * 1:length(transformation), lty = 2)
    
    # Plot 6
    plot(s, Y_hat[, i], type = "l", ylab = bquote(hat(Y)(t)))
    
  }
  
  return(cor(sunspots, Y_hat))
  
}

# Importing data and creating our dataset ---------------------------------

# Importing data
d <- read_csv("OrakeiBasin.csv")

# Adding Mo.inc/Mo.coh column
d <- d %>% mutate(Mo.ratio = Mo.inc / Mo.coh)

# Normalising all element columns
not_elements <- c(1:11, 62:75, 79, 88)
total_counts_elements <- rowSums(d[, -not_elements])
d <- d %>% mutate(total_count = total_counts_elements)
d[, -c(not_elements, 89)] <- d[, -c(not_elements, 89)] / (d[, "total_count"] %>% as.matrix() %>% as.numeric())

# Choosing start and end years (backwards from 2019)
start_year <- 1749
end_year <- 2000

# Creating sunspots time series
start_index <- 3177 - (2013 - start_year) * 12 + 9
end_index <- 3177 - (2013 - end_year) * 12 + 9
sunspots <- sunspot.month[start_index:end_index]

elements <- c(12, 13, 15:18, 20:38, 44:58, 61, 62, 64:72, 88)

# Creating mud time series
mm_per_year <- 0.677919

mm_needed <- (end_year - start_year) * mm_per_year
min_depth <- min(d$Depth) + (2016 - end_year) * mm_per_year / 1000
max_depth <- min_depth + mm_needed / 1000

element_data <- d %>% filter(Depth >= min_depth & Depth <= max_depth) %>% select(elements)
mo.coh_data <- d %>% filter(Depth >= min_depth & Depth <= max_depth) %>% select(Mo.coh)
mo.inc_data <- d %>% filter(Depth >= min_depth & Depth <= max_depth) %>% select(Mo.inc)
data <- element_data / as.numeric(as.matrix((mo.coh_data + mo.inc_data)))
data <- data * max(sunspots, na.rm = TRUE) / max(data)
data <- as.matrix(data)

s <- 1:length(sunspots)
t <- d %>% filter(Depth >= min_depth & Depth <= max_depth) %>% select(Depth) %>% as.matrix() %>% as.numeric()


# Results -----------------------------------------------------------------

# Mud data

n_knots <- 4
par(cex.lab = 1.5, cex.axis = 1.5, mar = c(5, 5, 2, 2))
(result <- time_warp_mv(get_a(n_knots),
  s, t, sunspots, data,
  parscale = get_parscale_mv(seq(-0.1, 0.1, by = 0.001), n_knots = n_knots),
  lambda_1 = 1, # if lambda are changed, make sure to change it in parscale and re-source the code
  lambda_2 = 1,
  lambda_3 = 1,
  i = 2))
sum(abs(result)) / length(result)

