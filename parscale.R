########################################
#
# Investigating parscale option for optim 
# 17/10/2019
#
########################################

hs <- seq(-0.1, 0.1, by = 0.001)

get_a <- function(n_knots) {
  
  a_i <- (max(t) - min(t)) / (max(s) - min(s)) * length(s) / n_knots
  
  return(c(min(t), rep(a_i, n_knots)))
  
}

get_res_mv <- function(h, i) {
  
  a <- get_a(n_knots)
  
  a[i] <- a[i] + h
  
  spline <- phi(s, a, order) 
  t_hat <- spline$t_hat
  t_hat_matched <- match_depths(t_hat, t)
  transformation <- sapply(t_hat_matched, function(x) {which(t == x)})
  Y_hat <- data[transformation, ]
  t_warped <- t[transformation]
  
  difference_penalty <- sum(diff(a[-1])^2)
  decrease_penalty <- sum(transformation %>% diff() <= 0)
  
  # Penalty for major deviation from average mm per year
  slope_penalty <- sum((a / 12 - 0.677919 / 1000)^2)
  
  result <- (sum(abs(cor(sunspots, Y_hat)) - 1 * difference_penalty - 1 * slope_penalty - 1 * decrease_penalty)) * -1
  
  return(result)
  
}

get_coef_mv <- function(i) {
  
  print(i)
  result <- sapply(hs, get_res_mv, i = i)
  model <- lm(result ~ hs, data = data.frame(result, hs))
  coef(model)[2]
  
}

get_parscale_mv <- function(hs = hs, t = t, data = data, s = s, sunspots = sunspots, n_knots = 4) {
  
  result <- sapply(1:(n_knots + 1), get_coef_mv)
  print(result)
  return(result)
  
}
