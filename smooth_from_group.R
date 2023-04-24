
##------------- FUNCTIONS FOR GROUPED DATA TO SMOOTH DIST ----------------------
## 
##
##  Adapted code from Rizzi et al. (2015) AJE
##
## 
##  Author: Benjamin Schlüter
##  Date: April 2023
##
##  Notes:
##  1) Smoothed mx increase "linearly", upper age bound. No Kannisto shape.
##     --> Combine e_j and Kannisto with s_hat()?
##  2) Computation by States will involve super small counts and pclm() 
##     might fail (focus on B&W?)
##  3) Smooth rates take smoothed exposures as fixed
##  4) Big population sizes makes the pclm fail to converge
##     but reducing pop size increase the smoothing
##     --> Reduce lambda range for large population size (ie white)
##  5) mx.s cannot be >2 otherwise qx>1 -> condition should be added
##
## -----------------------------------------------------------------------------


## Unmodified code from authors ------------------------------------------------

# Demo of the penalized composite link model (PCLM) for grouped counts
pclm <- function(y, C, X, lambda = 1, deg = 2, show = F){
        # Fit a PCLM (estimate b in ) E(y) = C %*% exp(X %*% b)
        # y = the vector of observed counts of length i
        # C = the composition matrix of dimension IxJ
        # X = the identity matrix of dimension JxJ; or B-spline basis
        # lambda = smoothing parameter
        # deg = order of differences of the components of b
        # show = indicates whether iteration details should be shown
        # Fit the penalized composite link model
        # Some preparations
        nx <- dim(X)[2]
        D <- diff(diag(nx), diff=deg)
        la2 <- sqrt(lambda)
        it <- 0
        bstart <- log(sum(y) / nx);
        b <- rep(bstart, nx);
        
        # Perform the iterations
        for (it in 1:50) {
                b0 <- b
                eta <- X %*% b
                gam <- exp(eta)
                mu <- C %*% gam
                w <- c(1 / mu, rep(la2, nx - deg))
                Gam <- gam %*% rep(1, nx)
                Q <- C %*% (Gam * X)
                z <- c(y - mu + Q %*% b, rep(0, nx - deg))
                Fit <- lsfit(rbind(Q, D), z, wt = w, intercept = F)
                b <- Fit$coef
                db <- max(abs(b - b0))
                if (show) cat(it, " ", db, "\n")
                if (db < 1e-6) break
                }
        cat(it, " ", db, "\n")
        
        # Regression diagnostic
        R <- t(Q) %*% diag(c(1 / mu)) %*% Q
        H <- solve(R + lambda * t(D) %*% D) %*% R
        fit <- list()
        fit$trace <- sum(diag(H))
        ok<-y>0&mu>0
        fit$dev <- 2 * sum(y[ok] * log(y[ok] / mu[ok]))
        fit$gamma <- gam
        fit$aic<- fit$dev + 2 * fit$trace
        fit$mu <- mu
        fit
}

## Example from paper:

# Simulate latent data
# m <- 130
# x <- 1:m
# xmean <- 80
# xsd <- 10
# set.seed(2012)
# f <- dnorm(x, xmean, xsd)
# mu <- 1e4 * f
# z <- rpois(m, lambda = mu)

# Compute the group counts
# gr <- seq(1, 86, by = 5)
# bnd <- 114
# ilo <- c(gr, bnd)
# ihi <- c(seq(5,85,5),115,130)
# n <- length(ihi)
# ihi[n] <- m
# y <- g <- 0 * ihi
# for (i in 1:n) {
#         y[i] <- sum(z[ilo[i]:ihi[i]])
#         }
# y <- c(y[1:17],sum(y[18:19]),0)
# for (i in 1:n) {
#         g[i] <- y[i] / (ihi[i] - ilo[i] + 1)
# }


# Make C matrix and (trivial) basis B
# C <- matrix(0, n, m)
# C[1:17, 1:85] <- kronecker(diag(17), matrix(1, 1, 5))
# C[18, 86:115] <- 1
# C[19, 116:130] <- 1
# 
# B <- diag(m)

# Solve PCLM
# lambda <- 10^7
# mod <- pclm(y, C, B,lambda = lambda, deg = 2)
# cat("lambda, ED & AIC:", lambda, mod$trace, mod$aic, "\n")

# # Plot data and fit
# plot(x, mu, type = "l", xlim = c(50, m), ylim = range(mod$gamma), lwd = 2,
#      xlab = "Age", ylab = "Number of events")
# lines(x, z, col = "darkgreen")
# points(x, z, pch=20,cex=0.8, col="darkgreen")
# lines(ihi, g, type = "S", col = "blue")
# lines(x, mod$gamma, col = "orange", lwd = 2)



## Define function using pclm() ------------------------------------------------

## This function works with group_by() and group_modify(),
## allowing to run it on all strata of interest.
## It uses the R code from Rizzi et al. (2015) AJE to extend
## mortality rates above age 85+.
## It returns a data set with smooth distribution above age 85.

## df needs to have age classes in df$age,
## population counts in df$pop, and
## death counts in df$n_deaths to work correctly.

## Argument for lambda could be added
smooth_from_group <- function(df, age.bnd = 110, rates = F) {
        
        
        ## 1ST SECTION: EXTRAPOLATE POPULATION COUNTS WITH PCML()
        
        ## y = the vector of population counts of length I.
        ## Complement histogram with zero death count
        ## following description in paper.
        y <- c(df$pop, 0)
        ## I consists of age classes in data + selected 
        ## upper bound (age.bnd) + interval set to zero
        age.class <- c(tail(df$age, -1), age.bnd, age.bnd + 10)
        I <- length(age.class)
        ## J consists of continuous age grid
        J <- length(0:max(age.class))
        
        ## Age starts at zero so row != Age (ie row=1 for Age=0)
        ## These indexes make the code more "readable"
        row.open.int <- max(df$age) + 1
        row.age.bnd <- age.bnd + 1
        
        ## C = the composition matrix of dimension IxJ
        ## If initial data has 5-year age classes,
        ## main modification for this function to work
        ## should be done in C.
        C <- matrix(0, I, J)
        ## 1-year age classes
        diag(C[1:(I-2), 1:(I-2)]) <- 1
        ## Extrapolate ages from 85 to age.bnd
        C[I-1, (I-1):(age.bnd+1)] <- 1
        ## Ages belonging to interval set to zero
        C[I, (age.bnd+2):J] <- 1
        
        ## X = the identity matrix of dimension JxJ
        B <- diag(J)
        
        ## Iterate on lambda to minimize AIC
        lambdas <- 10^seq(-5, 5, 0.1)
        nl <- length(lambdas)
        out <- list()
        for (i in 1:nl) {
                
                out[[i]] <- pclm(y, C, B,lambda = lambdas[i], deg = 2)
        }
        ## Get all AIC
        AIC <- unlist(lapply(out, "[[", "aic"))
        ## Get gamma (=extrapolated pop counts) associated with
        ## lambda minimizing AIC
        e_j <- c(out[[which.min(AIC)]]$gamma)
        ## Remove the ages where we set the pop counts to zero
        ## and keep extrapolated population counts
        pop.ext <- e_j[row.open.int:row.age.bnd]
        ## Make sure the sum of extrapolated pop and Inf_N_85 are equal
        cst.mult <- df$pop[row.open.int]/sum(pop.ext)
        pop.ext <- cst.mult * pop.ext
        
        
        
        ## 2ND SECTION: EXTRAPOLATE MX WITH PCML() AND GET DEATH COUNTS
        
        ## y = the vector of observed death counts of length i.
        y.rate <- df$n_deaths
        
        ## Each column j is multiplied by the respective e_j
        ## (obtained pop.ext)
        ## such that mu = C*gamma can be written.
        ## Gamma thus consist of rates.
        ## Note that uncertainty in e_j is not accounted for in 
        ## this stage. 
        C.rate <- C %*% diag(e_j)
        
        ## Modify C as death vector stops at open age interval
        ## (We don't add a hist. bin equal to zero as in pop case)
        C.rate <- C.rate[1:(I-1), 1:(age.bnd + 1)]
        
        ## Keep B part that matters
        # B.rate <- B[1:row.age.bnd, 1:row.age.bnd]
        ## Here we use B-spline instead to smooth mx for smaller pop
        ## These knots were obtained by trial and error (need both flexibility 
        ## for whites but robustness for ethnicity being less populous)
        B.rate   <- splines::bs(1:row.age.bnd, knots=c(0,1,5,15,30,40,70), degree=3)
        
        ## Iterate on lambda to minimize AIC
        ## Increase range with B-splines
        lambdas <- 10^seq(-7, 7, 0.1)
        nl <- length(lambdas)
        out <- list()
        for (i in 1:nl) {
                
                out[[i]] <- pclm(y.rate, C.rate, B.rate,lambda = lambdas[i], deg = 2)
        }
        ## Get all AIC
        AIC <- unlist(lapply(out, "[[", "aic"))
        ## Get gamma (=rates with extrapolation) associated with
        ## lambda minimizing AIC
        mx.pclm <- c(out[[which.min(AIC)]]$gamma)
        ## Death counts using extrapolated mx
        d_j <- mx.pclm * e_j[1:row.age.bnd]
        ## Keep extrapolated death counts
        dth.ext <- d_j[row.open.int:row.age.bnd]
        ## Make sure the sum of extrapolated death and Inf_D_85 are equal
        cst.mult <- df$n_deaths[row.open.int]/sum(dth.ext)
        dth.ext <- cst.mult * dth.ext
        
        ## Store outputs
        df.out <- tibble(age = 0:age.bnd,
                         pop.ext = c(head(df$pop, - 1), pop.ext),
                         dth.ext = c(head(df$n_deaths, - 1), dth.ext),
                         dth.pclm = d_j[1:row.age.bnd],
                         mx = c(df$mx, rep(NA, age.bnd - max(df$age))),
                         mx.pclm = mx.pclm[1:row.age.bnd])
        
        return(df.out)
}
