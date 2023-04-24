
##------------- FUNCTIONS TO SPLIT DEATH IN OPEN AGE INTERVAL ------------------
## 
##
##  Methods protocol
##  for the Human Mortality Database (V6).
##
##  The data consists of death and population counts by 1-year age class.
##  However, last open age class is 85+. The aim of the functions below is to
##  extrapolate death counts above age 84 yo by 1-year age classes.
##  
## 
##  Author: Benjamin Schlüter
##  Date: April 2023
##
##  Notes:
##  1) 
##
## -----------------------------------------------------------------------------


## Functions associated with Kannisto model ------------------------------------

S_kannisto <- function(par, x) {
        
        S_x <- ( (1 + par[1]) / (1 + (par[1]*exp(par[2]*x))) )^(1/par[2])
        return(S_x)
}

mu_kannisto <- function(par, x) {
        
        mu_x <- (par[1] * exp(par[2]*x)) / (1 + (par[1] * exp(par[2]*x)))
        return(mu_x)
}


## Function extrapol_kannisto ---------------------------------------------

## This function is defined to work with group_by() and group_modify(),
## allowing to run it on all strata of interest (i.e sex*year*race).

## df needs to have age classes in df$age,
## and population counts in df$pop,
## by 1-year age class.

## Argument for lambda could be added

extrapol_death_kannisto <- function(df, age.bnd = 100) {
        
        ## First row is 0 yo
        open.int <- max(df$age) + 1
        dth.above.65 <- sum(df$n_deaths[(open.int - 20):open.int])
        ## Fictitious survival fct
        ## S=1 by definition at age 65
        S <- c()
        for (x in 0:20) {
                dth.above.xi <- sum(df$n_deaths[(open.int - x):open.int]) 
                S[x+1] <- dth.above.xi/dth.above.65
        }
        S <- rev(S)
        
        ## Create data used in optim() to minimize residual sum of squares
        df.rss <- tibble(logS = tail(log(S), -1),
                         x = 1:20)
        ## Function minimizing the squared differences between the log
        ## of observed and predicted cumulative proportions.
        min.RSS <- function(data, par) {
                with(data, sum((( (1/par[2])*(log(1+par[1]) - log(1+(par[1]*exp(par[2]*x)))) ) - logS)^2))
        }
        ## Obtain initial values from a logistic regression (Feehan (2018)).
        ## Does not use mx at age 85+ 
        df.init <- tibble(mx.above65 = df$mx[(open.int - 20):(open.int-1)],
                          x = 0:19)
        par.init <- glm(mx.above65 ~ x, data = df.init, family = "binomial")$coef
        
        # df.init <- df.init %>%
        #         mutate( mx_hat = 1/(exp(-(par.init[1] + par.init[2]*x)) + 1) )
        # plot(x = df.init$x, y= df.init$mx.above65, t="p")
        # lines(x=df.init$x, y=df.init$mx_hat)
        
        ## Perform optimization and extract estimated parameters
        pars <- optim(par = c(exp(par.init[1]), par.init[2]), 
                      fn = min.RSS, 
                      data = df.rss,
                      method = "L-BFGS-B",
                      ## a and b have to be >0
                      lower = c(0, 0))$par
        
        # S_hat <- S_kannisto(pars, 0:20)
        # plot(x=0:20, y = S, t="p")
        # lines(x=0:20, y = S_hat, t="l")
        # mu_hat <- mu_kannisto(pars, 0:100)
        # plot(x=0:100, y = mu_hat, t="l")
        
        ## Predict d(x) and S(x)
        d_x <- c()
        S_x <- c()
        j <- 0
        for (x in seq(max(df$age), age.bnd+0.5, 0.5) - 65) {
                j <- j + 1
                d_x[j] <- S_kannisto(pars, x) - S_kannisto(pars, (x + 0.5))
                S_x[j] <- S_kannisto(pars, (x + 0.5))
        }
        ## Distribution of deaths by Lexis triangle within the open age interval 
        ## beginning at age 85
        D_l <- df$n_deaths[open.int] * (d_x[seq(1,31,2)]/S_kannisto(pars, 20))
        D_u <- df$n_deaths[open.int] * (d_x[seq(2,32,2)]/S_kannisto(pars, 20))
        
        ## Adjust proportionally (i.e., multiplied by a constant) so that their
        ## sum is equal to Inf_D_85.
        mult.cst <- df$n_deaths[open.int]/sum(D_l + D_u) 
        D.ext <- mult.cst*(D_l + D_u) 
        
        ## Distribution of exposure by Lexis triangle within the open age interval 
        ## beginning at age 85
        E_l <- df$pop[open.int] * (S_x[seq(1,31,2)]/S_kannisto(pars, 20))
        E_u <- df$pop[open.int] * (S_x[seq(2,32,2)]/S_kannisto(pars, 20))
        
        ## Adjust proportionally (i.e., multiplied by a constant) so that their
        ## sum is equal to Inf_E_85.
        mult.cst <- df$pop[open.int]/sum(E_l + E_u) 
        E.ext <- mult.cst * (E_l + E_u) 
        
        ## Store in df
        pop.out.K <- c(df$pop[1:max(df$age)], E.ext)
        dth.out.K <- c(df$n_deaths[1:max(df$age)], D.ext)
        mx.out.K <- dth.out.K/pop.out.K
        
        df.out <- tibble(age = 0:age.bnd,
                         pop.K = pop.out.K,
                         dth.K = dth.out.K,
                         mx.K = mx.out.K) 
        
        return(df.out)
        
        ## Check outputs (useful at function building stage/debugging)
        
        # df.out %>%
        #         ggplot(aes(x = age, y = pop.K)) +
        #         geom_line() +
        #         theme_bw()
        # df.out %>%
        #         ggplot(aes(x = age, y = dth.K)) +
        #         geom_line() +
        #         geom_point(aes(y = ifelse(age<=84, dth.K, NA))) +
        #         theme_bw()
        # df.out %>%
        #         ggplot(aes(x = age, y = mx.K)) +
        #         geom_line() +
        #         theme_bw() +
        #         scale_y_log10()
}

