library(data.table)
library(ineq)
library(Matrix)
library(recosystem)
library(tidyverse)
set.seed(305692189)


# Importing and cleaning data


options(timeout = 600)
url <- "http://konect.cc/files/download.tsv.libimseti.tar.bz2"
download.file(url, destfile = "libimseti.tar.bz2", mode = "wb")
untar("libimseti.tar.bz2", exdir = "libimseti_data")

list.files("libimseti_data", recursive = TRUE)

ratings_dt <- fread(
  "libimseti_data/libimseti/out.libimseti",
  skip = 1,
  header = FALSE,
  col.names = c("from_user", "to_user", "rating")
)

given_cnt <- ratings_dt[, .(n_given = .N), by = .(user_id = from_user)]
recv_cnt  <- ratings_dt[, .(n_recv  = .N), by = .(user_id = to_user)]
act <- merge(given_cnt, recv_cnt, by = "user_id")[n_given >= 20 & n_recv >= 20]
act[, total := n_given + n_recv]

top_ids <- act[order(-total)][1:2000, user_id]
sub_r   <- ratings_dt[from_user %in% top_ids & to_user %in% top_ids]

u_idx <- match(sub_r$from_user, top_ids)
v_idx <- match(sub_r$to_user, top_ids)
A <- sparseMatrix(i = u_idx, j = v_idx, x = 1, dims = c(2000, 2000))
A_sym <- A + t(A)

ev <- eigen(as.matrix(A_sym), symmetric = TRUE)
v_min <- ev$vectors[, 2000]

side_X <- top_ids[order(v_min)][1:1000]
side_Y <- top_ids[order(-v_min)][1:1000]

hetero_ratings <- sub_r[
  (from_user %in% side_X & to_user %in% side_Y) |
    (from_user %in% side_Y & to_user %in% side_X)
]

cat(sprintf(
  "Cross-side interaction share: %.1f%% (%d cross-side ratings out of %d)\n",
  100 * nrow(hetero_ratings) / nrow(sub_r), nrow(hetero_ratings), nrow(sub_r)
))

hetero_ratings[, is_test := runif(.N) < 0.20, by = from_user]
train_dt <- hetero_ratings[is_test == FALSE]
test_dt  <- hetero_ratings[is_test == TRUE]



# Estimating Preferences

all_sub_users <- c(side_X, side_Y)
train_dt[, u_idx := match(from_user, all_sub_users)]
train_dt[, i_idx := match(to_user, all_sub_users)]

r <- Reco()
train_data <- data_memory(
  user_index = train_dt$u_idx,
  item_index = train_dt$i_idx,
  rating     = train_dt$rating,
  index1     = TRUE
)

r$train(
  train_data,
  opts = list(dim = 20, lrate = 0.05, costp_l2 = 0.05, costq_l2 = 0.05, niter = 30, verbose = FALSE)
)

grid_xy <- CJ(x_idx = 1:1000, y_idx = 1:1000)
pred_xy <- data_memory(user_index = grid_xy$x_idx, item_index = grid_xy$y_idx + 1000L, index1 = TRUE)
grid_xy[, r_hat := r$predict(pred_xy, out_memory())]
grid_xy[, p_xy  := plogis((r_hat - 6.5) / 1.5)]
P_xy <- matrix(grid_xy$p_xy, nrow = 1000, ncol = 1000, byrow = TRUE)

grid_yx <- CJ(y_idx = 1:1000, x_idx = 1:1000)
pred_yx <- data_memory(user_index = grid_yx$y_idx + 1000L, item_index = grid_yx$x_idx, index1 = TRUE)
grid_yx[, r_hat := r$predict(pred_yx, out_memory())]
grid_yx[, p_yx  := plogis((r_hat - 6.5) / 1.5)]
P_yx <- matrix(grid_yx$p_yx, nrow = 1000, ncol = 1000, byrow = TRUE)



# Computing Reciprocal Fusion Scores

P_yx_T <- t(P_yx)

phi_arith <- (P_xy + P_yx_T) / 2
phi_geom  <- sqrt(P_xy * P_yx_T)
phi_harm  <- (2 * P_xy * P_yx_T) / (P_xy + P_yx_T + 1e-12)

# Implementing IPFP Solver for TU Matching Equilibrium

run_ipfp <- function(P_xy, P_yx_T, beta = 5, tol = 1e-6, max_iter = 1000) {
  P_tilde <- exp(beta * (P_xy + P_yx_T) / 2)
  n_x <- nrow(P_tilde)
  n_y <- ncol(P_tilde)
  
  s_x <- rep(1, n_x)
  s_y <- rep(1, n_y)
  
  for (iter in 1:max_iter) {
    s_x_old <- s_x
    
    a_x <- 0.5 * as.vector(P_tilde %*% s_y)
    s_x <- sqrt(1 + a_x^2) - a_x
    
    b_y <- 0.5 * as.vector(t(P_tilde) %*% s_x)
    s_y <- sqrt(1 + b_y^2) - b_y
    
    if (max(abs(s_x - s_x_old)) < tol) {
      cat(sprintf("IPFP (beta = %.1f) converged in %d iterations.\n", beta, iter))
      break
    }
  }
  
  P_tilde * outer(s_x, s_y)
}

mu_beta1 <- run_ipfp(P_xy, P_yx_T, beta = 1)
mu_beta5 <- run_ipfp(P_xy, P_yx_T, beta = 5)



# Semi-Synthetic Top-K Replay Simulation & Evaluation


rnoise <- function(n, dist = c("gumbel", "normal", "t"), df = 3) {
  dist <- match.arg(dist)
  switch(
    dist,
    gumbel = -log(-log(runif(n))),
    normal = rnorm(n, mean = 0, sd = pi / sqrt(6)),
    t      = rt(n, df = df)
  )
}

obs_like_rate <- mean(test_dt$rating >= 7)

### Simulating latent realized utilities u_{x,y} = logit(p_{x,y}) + eps_{x,y}
u_xy <- qlogis(pmin(pmax(P_xy, 1e-6), 1 - 1e-6)) + rnoise(1000 * 1000, dist = "gumbel")
u_yx <- qlogis(pmin(pmax(P_yx_T, 1e-6), 1 - 1e-6)) + rnoise(1000 * 1000, dist = "gumbel")

L_xy <- (u_xy > quantile(u_xy, 1 - obs_like_rate)) * 1L
L_yx <- (u_yx > quantile(u_yx, 1 - obs_like_rate)) * 1L

### Overwriting simulated pairs with ground-truth held-out test ratings (rating >= 7 => Like)
test_xy <- test_dt[from_user %in% side_X & to_user %in% side_Y]
test_xy[, `:=`(xi = match(from_user, side_X), yi = match(to_user, side_Y), like = as.integer(rating >= 7))]
L_xy[cbind(test_xy$xi, test_xy$yi)] <- test_xy$like

test_yx <- test_dt[from_user %in% side_Y & to_user %in% side_X]
test_yx[, `:=`(xi = match(to_user, side_X), yi = match(from_user, side_Y), like = as.integer(rating >= 7))]
L_yx[cbind(test_yx$xi, test_yx$yi)] <- test_yx$like

### Evaluating across K = 10, 20, 50
methods <- list(
  "Arithmetic Mean"        = phi_arith,
  "Geometric Mean"         = phi_geom,
  "Harmonic Mean"          = phi_harm,
  "TU Matching (beta = 1)" = mu_beta1,
  "TU Matching (beta = 5)" = mu_beta5
)

eval_topk <- function(score_mat, K) {
  rec_X <- matrix(0L, 1000, 1000)
  for (i in 1:1000) rec_X[i, order(-score_mat[i, ])[1:K]] <- 1L
  
  rec_Y <- matrix(0L, 1000, 1000)
  for (j in 1:1000) rec_Y[order(-score_mat[, j])[1:K], j] <- 1L
  
  shown_likes_X  <- rec_X * L_xy
  shown_likes_Y  <- rec_Y * L_yx
  mutual_matches <- shown_likes_X * shown_likes_Y
  
  incoming_likes <- c(colSums(shown_likes_X), rowSums(shown_likes_Y))
  user_matches   <- c(rowSums(mutual_matches), colSums(mutual_matches))
  
  list(
    summary = data.table(
      K              = K,
      Mutual_Matches = sum(mutual_matches),
      Gini_Likes     = Gini(incoming_likes),
      Coverage       = mean(user_matches > 0),
      Unmatched_Rate = mean(user_matches == 0)
    ),
    incoming_likes = incoming_likes
  )
}

K_grid <- c(10, 20, 50, 100)
eval_outputs <- lapply(names(methods), function(m_name) {
  lapply(K_grid, function(k) {
    out <- eval_topk(methods[[m_name]], k)
    out$summary[, Method := m_name]
    out$dist <- data.table(Method = m_name, K = k, Incoming_Likes = out$incoming_likes)
    out
  })
})

results_table <- rbindlist(lapply(eval_outputs, function(x) rbindlist(lapply(x, `[[`, "summary"))))
likes_dist_dt <- rbindlist(lapply(eval_outputs, function(x) rbindlist(lapply(x, `[[`, "dist"))))

print(results_table[, .(
  Method,
  K,
  Mutual_Matches,
  Gini_Likes     = round(Gini_Likes, 3),
  Coverage       = round(Coverage, 3),
  Unmatched_Rate = round(Unmatched_Rate, 3)
)])


# Data Visualization

ggplot(results_table, aes(x = Gini_Likes, y = Mutual_Matches, color = Method, shape = factor(K))) +
  geom_point(size = 4) +
  geom_line(aes(group = Method), linewidth = 0.8) +
  labs(
    title    = "TU Matching vs. Standard Fusion Rules on LibimSeTi",
    subtitle = "Superstar Congestion and Total Mutual Matches across Top-K",
    x        = "Incoming Likes Gini Coefficient",
    y        = "Total Mutual Matches",
    shape    = "Top-K Limit"
  ) +
  theme_minimal(base_size = 12)

lorenz_dt_20 <- likes_dist_dt[K == 20, {
  x_sorted <- sort(Incoming_Likes)
  .(
    cumu_users = seq_along(x_sorted) / .N,
    cumu_likes = cumsum(x_sorted) / sum(x_sorted)
  )
}, by = Method]

ggplot(lorenz_dt_20, aes(x = cumu_users, y = cumu_likes, color = Method)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Lorenz Curve of Incoming Likes on LibimSeTi (Top-K = 20)",
    x        = "Cumulative Share of Users",
    y        = "Cumulative Share of Incoming Likes"
  ) +
  theme_minimal(base_size = 12)

lorenz_dt_100 <- likes_dist_dt[K == 100, {
  x_sorted <- sort(Incoming_Likes)
  .(
    cumu_users = seq_along(x_sorted) / .N,
    cumu_likes = cumsum(x_sorted) / sum(x_sorted)
  )
}, by = Method]

ggplot(lorenz_dt_100, aes(x = cumu_users, y = cumu_likes, color = Method)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Lorenz Curve of Incoming Likes on LibimSeTi (Top-K = 100)",
    x        = "Cumulative Share of Users",
    y        = "Cumulative Share of Incoming Likes"
  ) +
  theme_minimal(base_size = 12)
