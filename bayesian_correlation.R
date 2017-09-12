library(rstan)

load("processed_data_for_input.RData")

genes <- as.matrix(dat$genes)
volumes <- as.matrix(dat$volumes)

volumes_cor <- cor(volumes)

volumes_cor_lt <- volumes_cor[upper.tri(volumes_cor)]
volumes_cor_xf <- log ( (1 + volumes_cor_lt) / (1 - volumes_cor_lt) )

G <- nrow(genes)
P <- ncol(genes)
T <- length(volumes_cor_xf)

model <- stan_model("./bayesian_correlation_model.stan")

data <-
  list(genes = genes
     , volumes_cor_xf
     , G = G
     , P = P
     , T = T) 

stan_optim <- optimizing(model, data = data)
