source('function_h2.R')
n <- 2e5
h2.common <- 0.4
#h2.rare <- 0.2
#h2.env <- 1 - h2.rare - h2.common
prs1.e <- 0.1
prs2.e <- 10
beta.real <- 3
beta.assumed <- 3

symm <- TRUE
iter <- 1
h2.est1 <- matrix( nrow=iter, ncol=4 )

for( i in 1:iter ){
    common.effects <- sqrt(h2.common) * rnorm( n=n )
    rare.effects1 <- beta.real * rbinom(n=n,p=0.005,size=1) # captured rare effects
    rare.effects2 <- beta.real * rbinom(n=n,p=0.01,size=1) # missed rare effects
    if( symm ){
        rare.effects1 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects1 # captured rare effects
        rare.effects2 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects2 # missed rare effects
    }
    h2.env <- 1 - h2.common - var(rare.effects1) - var(rare.effects2)
    env <- sqrt(h2.env) * rnorm( n=n )

    y <- common.effects + rare.effects1 + rare.effects2 + env
    yy <- quantile.normalise(y)
    yyy <- yy

    prs1 <- common.effects + sqrt(prs1.e) * rnorm( n=n ) # less noise
    prs2 <- common.effects + sqrt(prs2.e) * rnorm( n=n ) # more noise

    prs1 <- prs1/sd(prs1)
    prs2 <- prs2/sd(prs2)

    ptr <- which( rare.effects1==0 & rare.effects2==0 )
    y.prime <- yyy[ptr]
    mu.y = mean(y.prime)
    sigma.y = sd(y.prime)
#    m.yprime <- 0
    fit <- summary(lm( yyy ~ I(prs1) ))
    test <- prs.test( (prs1), yyy )
    p.in.tail1 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex1 <- h2.rare.big( p.in.tail1, beta=beta.assumed )

    fit <- summary(lm( yyy ~ I(prs2) ))
    test <- prs.test( (prs2), yyy )
    p.in.tail2 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex2 <- h2.rare.big( p.in.tail2, beta=beta.assumed )

    fit <- summary(lm( yyy ~ I(prs1+rare.effects1) ))
    test <- prs.test( (prs1+rare.effects1), yyy )
    p.in.tail3 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex3 <- h2.rare.big( p.in.tail3, beta=beta.assumed )

    fit <- summary(lm( yyy ~ I(prs2+rare.effects1) ))
    test <- prs.test( (prs2+rare.effects1), yyy )
    p.in.tail4 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex4 <- h2.rare.big( p.in.tail4, beta=beta.assumed )

#    print( c( p.in.tail1, p.in.tail2, p.in.tail3, p.in.tail4 ) )
    print( c( ex1$h2, ex2$h2 ) )
    h2.est1[i,] <- c( ex1$h2, ex2$h2 )
}
print( c( var(rare.effects1), var(rare.effects2), var(rare.effects1)+var(rare.effects2) ) )

fit <- summary(lm( yyy ~ I(prs1) ))
r <- sqrt(fit$r.squared)
ptr <- which( rare.effects1==0 & rare.effects2==0 )
#ptr <- 1:n
y.prime <- yyy[ptr]
prs.prime <- prs1[ptr]

fit <- summary( lm(prs.prime ~ y.prime ))
#print(fit)
#print( sqrt(fit$r.squared) / sd(y.prime) )
#print( sqrt(fit$r.squared) * mean(y.prime) / sd(y.prime) )

K = 1
mu.y = mean(y.prime)
sigma.y = sd(y.prime)
K.prime = (K - mu.y)/sigma.y
#print( mean( prs.prime[ y.prime>K ] ) )
#print( sqrt(fit$r.squared) * mean( y.prime[ y.prime>K ] ) / sigma.y - mu.y / sigma.y )
#print( sqrt(fit$r.squared) * (mu.y + sigma.y*dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE)) / sigma.y - mu.y / sigma.y )
mu.y = 0
K.prime = (K - mu.y)/sigma.y
#print( sqrt(fit$r.squared) * (mu.y + sigma.y*dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE)) / sigma.y - mu.y / sigma.y )
mu.y = mean(y.prime)
sigma.y = 1
K.prime = (K - mu.y)/sigma.y
#print( sqrt(fit$r.squared) * (mu.y + sigma.y*dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE)) / sigma.y - mu.y / sigma.y )
mu.y = 0
sigma.y = 1
K.prime = (K - mu.y)/sigma.y
#print( sqrt(fit$r.squared) * (mu.y + sigma.y*dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE)) / sigma.y - mu.y / sigma.y )

#print( r * ( mean( y.prime[y.prime>2] ) - mean(y.prime) ) )
#print( mean( prs.prime[y.prime>2] ) )

test <- prs.test( prs1, yy )
prs.test( prs2, yy )
prs.test( (prs2+rare.effects1), yy )

