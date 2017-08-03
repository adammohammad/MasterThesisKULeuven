# Master Thesis KU Leuven
# Adam Mohammad
# Title: A Study on the Effect of Misspecifying the Association Structure in Dynamic Predictions 
# Obtained from Joint Models for Longitudinal and Survival Data
# August 2017
# Code for simulating data from the joint model with the current value and slope association structure

#### SIMULATE DATA ####
n=N
# at which time points longitudinal measurements are supposed to be taken
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t_max))))) 
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
covariate <-runif(n,-2,2)
DF <- data.frame(year = times, group = factor(rep(group, each = K)),covariate=rep(covariate, each=K))
# design matrices for the longitudinal measurement model
X <- model.matrix(~ 0 + group +covariate+ group:ns(year, knots = Int_kn, Boundary.knots = Bound_kn), data = DF)
Z <- model.matrix(~ ns(year, knots = Int_kn, Boundary.knots = Bound_kn), data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group,"Covariate"=covariate)

#simulate random effects#
b <- mvrnorm(n, rep(0, nrow(D)), D)

# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
y <- rnorm(n * K, eta.y, sigma)

# simulate event times
eta.t <- as.vector(W %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group0 <- 1 - group[i]
    group1 <- group[i]
    Covariate <- covariate[i]
    BS <- ns(s, knots = Int_kn,  Boundary.knots = Bound_kn)
    XX <- cbind(group0, group1, Covariate, group0*BS[, 1], group1*BS[, 1],
                group0*BS[, 2], group1*BS[, 2], group0*BS[, 3], group1*BS[, 3])
    ZZ <- cbind(1, BS)
    f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
    ###
    dBS <- dns(s, knots = Int_kn, Boundary.knots = Bound_kn)
    XXd <- cbind(group0*dBS[, 1], group1*dBS[, 1],
                 group0*dBS[, 2], group1*dBS[, 2], group0*dBS[, 3], group1*dBS[, 3])
    ZZd <- dBS
    f2 <- as.vector(XXd %*% betas[4:9] + rowSums(ZZd * b[rep(i, nrow(ZZd)), 2:4]))
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha1 + f2 * alpha2)
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
  Up <- 50
  tries <- 5
  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up + 200
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}

na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
group <- group[na.ind]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from a Beta distribution, with the mean depending on group
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- numeric(n)
Ctimes[group == 0] <- runif(sum(group == 0), 0, 2 * meanCens0)
Ctimes[group == 1] <- runif(sum(group == 1), 0, 2 * meanCens1)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################

# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
names(dat) <- c("time", "group","Covariate", "id", "y", "TimeEv", "event")
summary(Time)

set <- sample(unique(id), 500)
train_data <- dat[!dat$id %in% set, ]
train_data$id <- match(train_data$id, unique(train_data$id))
test_data <- dat[dat$id %in% set, ]
test_data$id <- match(test_data$id, unique(test_data$id))


# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W,
   eta.t, eta.y, trueTimes, u, Root, invS, set, dat, times, group, i, tries, Up,  DF)