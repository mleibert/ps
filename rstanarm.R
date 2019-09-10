
training <- Loans[1:10000,]
testing <- Loans[-c(1:10000),]

with(training,table(term,home_ownership,emp_length != 0))


library(rstanarm)
options(mc.cores=  parallel::detectCores())



post <- stan_glm(y ~ log(loan_amnt) + term + home_ownership +
                   log(annual_inc) + I(emp_length != 0), 
                 family = binomial(link = "logit"), data = training,
                 prior = normal(), prior_intercept = normal(), QR = TRUE)
post

head(bh)
bh

post <- stan_glm( cmedv~. , 
                  family = gaussian, data = bh,
                  prior = laplace(0,.1),   QR = TRUE)

post3 <- stan_glm( cmedv~. , 
                 family = gaussian, data = bh,
                 prior = laplace(1,.1),   QR = TRUE)

post4 <- stan_glm( cmedv~. , 
                  family = gaussian, data = bh,
                  prior = laplace(2,.1),   QR = TRUE)

post2 <- stan_glm( cmedv~. , 
                  family = gaussian, data = bh,
                      QR = TRUE)


df1 <- as.matrix(post)
df2 <- as.matrix(post2)
df3 <- as.matrix(post3)
df4 <- as.matrix(post4)

colno = "age"
colno = which(colnames(bh) == colno)

Z <- data.frame(rr=df1[,colno],unpenalized=df2[,colno])
dendata <- melt(Z)
ggplot(dendata, aes(x=value, fill=variable)) + geom_density(alpha=0.25)

Z <- data.frame(rr=df3[,colno],unpenalized=df2[,colno])
dendata <- melt(Z)
ggplot(dendata, aes(x=value, fill=variable)) + geom_density(alpha=0.25)

Z <- data.frame(lasso1=df1[,colno],unpenalized=df2[,colno],lasso2=df3[,colno],lasso3=df4[,colno])
dendata <- melt(Z)
ggplot(dendata, aes(x=value, fill=variable)) + geom_density(alpha=0.25)
 

#### #### #### #### #### 

B <- 5000
X <- model.matrix(cmedv~. , data = bh)
y <- bh$cmedv
p <- ncol(bh)
n = nrow(bh)

#### 


betamat <- matrix(NA, B, p )
betamat[1,] <-  coef(fit)
sigs <- rep(NA,B)
sigs[1] <- sigma(fit)^2


for (i in 2:B) { 
 
  e = y-X%*%betamat[i-1,]
  sigs[i] <- rinvgamma(1, n/2,  .5 * t(e) %*% e )
  
  betamat[i,] <-  rmvnorm(1, coef(fit) , sigs[i] * xtxi ) 
  
  # Sample Sigma^2 #
  sigs[i] <- rinvgamma(1, n/2, sum((y-X%*%betamat[i,])^2)/2)
}


betamat <- tail(betamat,B-1000)
colnames(betamat) <- colnames(X)
tail(betamat)

#### #### #### #### 

Betas <- taus <- matrix(NA, B, p )
ss <- sigs <- rep(NA,B)

fit <- lm(cmedv~. , data = bh)

Lambda = 0.3

Betas[1,] <- betamat[1,] <-  coef(fit)
ss[1] <- sigs[1] <- sigma(fit)^2
taus[1,] <- 1/rinvgauss( p ,  sqrt(Lambda*ss[1]/Betas[1,]^2) , Lambda^2 )

xtx <- t(X) %*% X
xtxi <- solve(xtx)
txy <- t(X) %*% y
ShapeParameterIG <- (n-1)/2 + p/2 + 1


for( i in 2:B){
  
  Dt <- solve( diag(taus[i-1,]) )
  A <- solve( xtx +  Dt  )
  Betas[i,] <- rmvnorm(1 , A %*% txy, ss[i-1] * A )
  ss[i] <- 1/rgamma( 1, ShapeParameterIG , 
                     .5 * t(y - X %*% Betas[i,] ) %*%  (y - X %*% Betas[i,] ) +
                       .5 * t(Betas[i,]) %*% Dt %*% Betas[i,] + 1 )
  
  taus[i,] <- 1/rinvgauss( p ,  sqrt(Lambda*ss[i]/Betas[i,]^2) , Lambda^2 )
  
}


 

Betas <- tail(Betas,B-1000)
colnames(Betas) <- colnames(X)
tail(Betas)

#### #### #### #### 

colno = "age"
colno = which(colnames(bh) == colno)


Z <- data.frame(Lasso=Betas[,colno], unpenalized=betamat[,colno] )
dendata <- melt(Z)
ggplot(dendata, aes(x=value, fill=variable)) + geom_density(alpha=0.25)


