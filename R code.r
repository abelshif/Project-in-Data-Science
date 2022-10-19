##Rcode for data generation, parameter estimation,and
#bias/variance calculation####
set.seed(16)
N  = 5000
x1 = sample(1:3, N, replace=TRUE)
x2 = runif(N, 0.05, 0.5)
X  = matrix(cbind(x1,x2),nrow=N,ncol=2)
X  = cbind ( matrix(1,N,1),X ) # Intercept

### calculation of ln_e_0t ####
#let y_0t =ln_e_0t
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
v_t_given_u_t=rnorm(N, mu_v_t_given_u_t, sigma_v_t_given_u_t)
gamma = 1
gamma = as.vector(c(1,1,1))
z = X%*%gamma + v_t_given_u_t
threshhold = log(4.5)    #log(4.5)is in millions of Swedish Kronor
w=0
y_eut = ifelse(z < threshhold , w, z) 
censored = ( y_eut == 0 )
ncensored = sum(censored) 
X[censored,]
k1 = ncol(X)
#let y_t = ln_e_t
### reparametrized parameters ####
alpha = 1/sigma_u_t
tau = beta/sigma_u_t
omega = 1/sigma_v_t*a
gamma = gamma
theta = rho/a

### the likelihood function ####
im.exp = function(zeta,y_0t,y_eut,N,X)
{
  k1 = ncol(X)
  alpha = zeta[1]^2
  tau = zeta[2:(k1+1)]
  omega = zeta[k1+2]^2
  gamma = zeta[(k1+3):(k1+k1+2)]
  theta = zeta[k1+k1+3]
  censored = ( y_eut == 0 )
  ncensored = sum(censored)

  loglik_y_t=N*(-1/2)*log(2*pi)+N*log(alpha)-1/2*sum(((alpha*y_0t)-(X %*%tau))^2)+(N-ncensored)*(-1/2)*log(2*pi) 
             +(N-ncensored)*log(omega)-1/2*sum((omega*(y_eut[!censored]-(X[!censored,]%*%gamma))-theta*((alpha*y_0t[!censored]) 
             -(X[!censored,]%*%tau)))^2)+sum(log(pnorm(-1*omega*(X[censored,]%*%gamma)-(theta*alpha*y_0t[censored]) 
             +theta*(X[censored,]%*% tau)))) 

  return(-loglik_y_t)
}
im.exp_grad<- function(zeta,y_0t,y_eut,N,X) 
{
  k1 = ncol(X)
  alpha = zeta[1]^2
  tau = zeta[2:(k1+1)]
  omega = zeta[k1+2]^2
  gamma = zeta[(k1+3):(k1+k1+2)]
  theta = zeta[k1+k1+3]
  censored = ( y_eut == 0 )
  ncensored = sum(censored)

  lamda_t=(dnorm(1*zeta[k1+2]^2*(X[censored,]%*%zeta[(k1+3):(k1+k1+2)])-(zeta[k1+k1+3]*zeta[1]^2*y_0t[censored])+zeta[k1+k1+3]*
          (X[censored,]%*%zeta[2:(k1+1)])))/(pnorm(-1*zeta[k1+2]^2*(X[censored,]%*%zeta[(k1+3):(k1+k1+2)])-
          (zeta[k1+k1+3]*zeta[1]^2* y_0t[censored])+zeta[k1+k1+3]*(X[censored,] %*% zeta[2:(k1+1)])))

  FSH = (zeta[1]^2*y_0t-(X%*%zeta[2:(k1+1)]))
  JMM = (zeta[k1+2]^2*(y_eut[!censored]-(X[!censored,]%*%zeta[(k1+3):(k1+k1+2)]))-zeta[k1+k1+3]*
        ((zeta[1]^2*y_0t[!censored])-(X[!censored,] %*% zeta[2:(k1+1)])))
  
  grad=c(N*(2/zeta[1])-sum(2*zeta[1]*y_0t*FSH)+sum(2*zeta[1]*zeta[k1+k1+3]*y_0t[!censored]*JMM)+
       sum(-2*lamda_t*zeta[1]*zeta[k1+k1+3]*y_0t[censored]),
       t(X)%*% FSH + -1*zeta[k1+k1+3]*(t(X[!censored,] )%*%JMM)+zeta[k1+k1+3]*(t(X[censored,])%*%lamda_t),
       (N-ncensored)*(2/zeta[k1+2])-sum(2*zeta[k1+2]*JMM*(y_eut[!censored]-(X[!censored,]%*%zeta[(k1+3):(k1+k1+2)])))-
       2*zeta[k1+2]*crossprod(lamda_t,(X[censored,]%*%zeta[(k1+3):(k1+k1+2)])),
       zeta[k1+2]^2*(t(X[!censored,])%*%JMM)-zeta[k1+2]^2*(t(X[censored,])%*%lamda_t),
       sum(JMM*(2*zeta[1]*y_0t[!censored]-(X[!censored,]%*%zeta[2:(k1+1)])))+
       sum(-1*lamda_t*(2*zeta[1]*y_0t[censored]-(X[censored,]%*%zeta[2:(k1+1)])))) # nolint

  return(-grad)
}

opt_grad=optim(c(alpha,tau,omega,gamma,theta),y_0t=y_0t,y_eut=y_eut,
               X=X,N=N,im.exp,im.exp_grad,method="BFGS",
               control=list(maxit=400))
opt_grad
alpha_hat = opt_grad$par[1]
tau_hat =   opt_grad$par[2:(k1+1)]
omega_hat = opt_grad$par[k1+2]
gamma_hat = opt_grad$par[(k1+3):(k1+k1+2)]
theta_hat = opt_grad$par[k1+k1+3]
cat("", "Parameter estimates", c(alpha_hat, tau_hat, omega_hat, gamma_hat, theta_hat), "", sep = "\n") # nolint
#print(c(alpha_hat, tau_hat, omega_hat, gamma_hat, theta_hat))

### real estimates ####
sigma_u_t_hat = 1/(alpha_hat)^2
beta_hat= tau_hat* sigma_u_t_hat
gamma_hat = gamma_hat
rho_hat = sqrt(theta_hat^2/(1+theta_hat^2))
sigma_v_t_hat = 1/((omega_hat)^2*sqrt(1-rho_hat^2))

cat("","Real estimates", c(beta_hat, gamma_hat, rho_hat, sigma_u_t_hat, sigma_v_t_hat), "", sep = "\n")
#print(c(beta_hat, gamma_hat, rho_hat, sigma_u_t_hat, sigma_v_t_hat))

### parameter estimates for 500 replicates ####
result = matrix(0,500,9)
x1 = sample(1:3, N, replace=TRUE)
x2 = runif(N, 0.05, 0.5)
X = matrix(cbind(x1,x2),nrow=N,ncol=2)
X = cbind (matrix(1,N,1),X )
for(i in 1:500)
{
  sigma_u_t = 2
  mu_u_t = 0
  u_t = rnorm(N, mu_u_t, sigma_u_t)
  beta = 1
  beta = as.vector(c(1,1,1))
  y_0t = X%*%beta + u_t
  sigma_v_t = 3.5
  rho = rho
  a = sqrt(1-rho^2)
  mu_v_t_given_u_t = (sigma_v_t/sigma_u_t)*rho*u_t
  sigma_v_t_given_u_t = a*sigma_v_t
  v_t_given_u_t=rnorm(N, mu_v_t_given_u_t, sigma_v_t_given_u_t)
  gamma = 1
  gamma = as.vector(c(1,1,1))
  z = X%*%gamma + v_t_given_u_t
  threshhold = log(4.5)   #log(4.5)is in millions of Swedish Kronor
  w=0
  y_eut = ifelse(z < threshhold , w, z) 
  censored = ( y_eut == 0 )
  ncensored = sum(censored) 
  X[censored,]
  alpha = 1/sigma_u_t
  tau = beta/sigma_u_t
  omega = 1/(sigma_v_t*a)
  gamma = gamma
  theta = rho/a
  opt_grad = optim(c(alpha,tau,omega,gamma,theta),y_0t=y_0t,
                 y_eut=y_eut,X=X,N=N,im.exp, im.exp_grad, 
                 method ="BFGS",control=list(maxit=400))
  opt_grad
  result[i,1:9] = opt_grad$par
}
k1 = ncol(X)
alpha_hat_i = result[,1]
tau_hat_i = result[,2:(k1+1)]
omega_hat_i = result[,k1+2]
gamma_hat_i = result[,(k1+3):(k1+k1+2)]
theta_hat_i = result[,k1+k1+3]
#print(c(alpha_hat_i, tau_hat_i, omega_hat_i, gamma_hat_i, theta_hat_i))

###  real estimtes for 500 replicates  ####
res = matrix(0,500,9)
for(i in 1:500)
{
  res[,k1+k1+2] = 1/(alpha_hat_i)^2
  res[,1:k1] = tau_hat_i* res[,k1+k1+2]
  res[,(k1+1):(k1+k1)] = gamma_hat_i
  res[,k1+k1+1] = sqrt(theta_hat_i^2/(1+ theta_hat_i^2))
  res[,k1+k1+3] = 1/((omega_hat_i)^2*sqrt(1-res[,k1+k1+1]^2))
  res
} 
beta_hat_i = res[,1:k1]
gamma_hat_i = res[,(k1+1):(k1+k1)]
rho_hat_i = res[,k1+k1+1]
sigma_u_t_hat_i = res[,k1+k1+2]
sigma_v_t_hat_i = res[,k1+k1+3]
print(c(beta_hat_i, gamma_hat_i, rho_hat_i, sigma_u_t_hat_i, sigma_v_t_hat_i))

### Testing bias of real estimates for 500 replictaes ####
colMeans(beta_hat_i)-beta
colMeans(gamma_hat_i)- gamma
mean(rho_hat_i)- rho
mean(sigma_u_t_hat_i)- sigma_u_t
mean(sigma_v_t_hat_i)- sigma_v_t

### Testing variance of real estimates for 500 replicates ####
apply(beta_hat_i, 2, var)
apply(gamma_hat_i, 2, var)
var(rho_hat_i)
var(sigma_u_t_hat_i)
var(sigma_v_t_hat_i)

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
cat("", "OLS estimates", c(alpha0, tau0, omega0, gamma0, theta0), "", sep = "\n")
#print(c(alpha0, tau0, omega0, gamma0, theta0))

