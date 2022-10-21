### OLS estimates used as initial values ####
set.seed(16)
N = 5000
x1 = sample(1:3, N, replace=TRUE)
x2 = runif(N, 0.05, 0.5)
X = matrix(cbind(x1,x2),nrow=N,ncol=2)
X = cbind ( matrix(1,N,1),X ) # Intercept

### calculation of ln_e_0t ####
#let y_0t = ln_e_0t
sigma_u_t = 2
mu_u_t = 0
u_t = rnorm(N, mu_u_t, sigma_u_t)
beta = 1
beta = as.vector(c(1,1,1))
y_0t = X%*%beta + u_t

### calculation of ln_e_eut ####
#let y_eut = ln_e_eut
sigma_v_t = 3.5
rho = 0.5
a = sqrt(1-rho^2)
mu_v_t_given_u_t = (sigma_v_t/sigma_u_t)*rho*u_t
sigma_v_t_given_u_t = a*sigma_v_t
v_t_given_u_t = rnorm(N, mu_v_t_given_u_t, sigma_v_t_given_u_t)
gamma = 1
gamma = as.vector(c(1,1,1))
z = X%*%gamma + v_t_given_u_t

# calculate the estimates BETA HAT
A = t(X)%*%X
Ainv = solve (A)
betahat = Ainv %*% ( t(X)%*%y_0t )
resl = y_0t - X%*%betahat
sse = crossprod( resl )
s_u = sqrt( sse/(N-3) )
betahat = as.vector(betahat)
s_u = as.vector(s_u)

# calculate the estimates GAMMA HAT
A = t(X)%*%X
Ainv = solve (A)
gammahat = Ainv %*% ( t(X)%*%z )
resl = z - X%*%gammahat
sse = crossprod( resl )
s_v = sqrt( sse/(N-3) )
gammahat = as.vector(gammahat)
s_v = as.vector(s_v)
resl1 = residuals(lm(y_0t ~ X))
resl2 = residuals(lm(z ~ X))
cor(resl1,resl2)
rhohat = cor(resl1,resl2)
rhohat
a = sqrt(1-rhohat^2)
a
alpha0 = 1/s_u
tau0   = betahat/s_u
omega0 = 1/s_v*a
gamma0 = gammahat
theta0 = rhohat/a

### initial values ####
cat("", "Initial values", sep = "\n")
print(c(alpha0, tau0, omega0, gamma0, theta0))
