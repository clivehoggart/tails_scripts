source('~/bin/functions.R')
source('functions.R')
source('function_h2.R')
n <- 2e5
n.sim <- 1e5
h2.common <- 0.1258
prs.r2.1 <- 0.025
prs.r2.2 <- 0.05
prs1.e <- h2.common * ( h2.common / prs.r2.1 - 1 )
prs2.e <- h2.common * ( h2.common / prs.r2.2 - 1 )

symm <- TRUE
iter <- 1
h2.est1 <- matrix( nrow=iter, ncol=6 )

rare.maf <- 1e-4

# Gives age of menopause signal
beta.real <- 4
beta.assumed <- 2
n.rare.effects1 <- 90
n.rare.effects2 <- 20

#beta.real <- 3
#beta.assumed <- 3
#n.rare.effects1 <- 50
#n.rare.effects2 <- 50

for( i in 1:iter ){
    common.effects <- sqrt(h2.common) * rnorm( n=n )
    rare.effects1 <- apply(matrix(ncol=n,data=replicate( n,
                                                        rnorm( n.rare.effects1, -beta.real, 0 ) *
                                                        rbinom(n=n.rare.effects1,p=rare.maf,size=2))), 2, sum )
    rare.effects2 <- apply(matrix(ncol=n,data=replicate( n,
                                                        rnorm( n.rare.effects2, beta.real, 0 ) *
                                                        rbinom(n=n.rare.effects2,p=rare.maf,size=2))), 2, sum )
#    if( symm ){
#        rare.effects1 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects1 # captured rare effects
#        rare.effects2 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects2 # missed rare effects
#    }

    h2.env <- 1 - h2.common - var(rare.effects1) - var(rare.effects2)
    env <- sqrt(h2.env) * rnorm( n=n )
    print(h2.env)

    y <- common.effects + rare.effects1 + rare.effects2 + env
    yy <- quantile.normalise(y)
    yyy <- yy

    prs1 <- common.effects + sqrt(prs1.e) * rnorm( n=n ) # less noise
    prs2 <- common.effects + sqrt(prs2.e) * rnorm( n=n ) # more noise
    delta.prs.r2 <- cor( prs1, yyy )^2 - prs.r2.1
    prs1.ee <- prs1.e / 10
    while( delta.prs.r2>0.0005 ){
        prs1 <- prs1 + sqrt(prs1.ee) * rnorm( n=n )
        delta.prs.r2 <- cor( prs1, yyy )^2 - prs.r2.1
    }
    prs1 <- prs1/sd(prs1)
    prs2 <- prs2/sd(prs2)

    ptr <- which( rare.effects1==0 & rare.effects2==0 )
    y.prime <- yyy[ptr]
    prs.prime <- prs1[ptr]
    mu.y = mean(y.prime)
    sigma.y = sd(y.prime)
#    m.yprime <- 0

    fit <- summary(lm( yyy ~ I(prs1) ))
    test <- prs.test( (prs1), yyy )

#    p.in.tail <- est.prop.in.tail.emp( test[,'effect'], beta=beta.assumed, r2=fit$r.squared, y.prime=y.prime, prs=prs.prime )
#    h2.rare.big.emp( p.in.tail, beta=beta.assumed, y.prime=y.prime, prs=prs.prime )
#    mean( yy>2.3 & rare.effects2>0 ) / mean( yy>2.3 )

    print(test)
    print(fit$r.squared)
    p.in.tail1 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex1 <- h2.rare.big2( p.in.tail1, beta=beta.assumed, rare.maf=1e-4, mu.y=mu.y, sigma.y=sigma.y )

    fit <- summary(lm( yyy ~ I(prs2) ))
    test <- prs.test( (prs2), yyy )
    p.in.tail2 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex2 <- h2.rare.big2( p.in.tail2, beta=beta.assumed, mu.y=mu.y, sigma.y=sigma.y )

    fit <- summary(lm( yyy ~ I(prs1) ))
    test <- prs.test( (prs1), yyy )
    ex3 <- h2.est.emp( n=n.sim, effect.size=test[,'effect'], prs.r2=fit$r.squared, h2.common=h2.common,
                      beta=beta.assumed, sd.beta=0 )

    fit <- summary(lm( yyy ~ I(prs2) ))
    test <- prs.test( (prs2), yyy )
    ex4 <- h2.est.emp( n=n.sim, effect.size=test[,'effect'], prs.r2=fit$r.squared, h2.common=h2.common,
                      beta=beta.assumed, sd.beta=0 )

#    print( c( p.in.tail1, p.in.tail2, p.in.tail3, p.in.tail4 ) )
    print( c( var(rare.effects1) + var(rare.effects2), 1-var(yy[ptr]), ex1$h2, ex2$h2 ) )
    print( apply(ex3,2,mean) )
    print( apply(ex4,2,mean) )
#    h2.est1[i,] <- c( var(rare.effects1) + var(rare.effects2), 1-var(yy[ptr]),
#                     ex1$h2, ex2$h2, apply(ex3,2,mean)['h2'], apply(ex4,2,mean)['h2'] )
}

#    mean( yyy>2.3 & rare.effects2>0 ) / ( mean( y.prime>2.3-3 ) * rare.maf * 40 )

#    fit <- summary(lm( yyy ~ I(prs1) ))
#    test <- prs.test( (prs1), yyy )
#    p.in.tail1 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
#    unlist(h2.rare.big2( p.in.tail1, beta=beta.assumed, rare.maf=1e-4, mu.y=mu.y, sigma.y=sigma.y ))

#    prs.prime <- prs1[ptr]
#    p.in.tail11 <- est.prop.in.tail.emp( test[,'effect'], beta.assumed, y.prime=y.prime, prs=prs.prime )
#    unlist(h2.rare.big.emp( p.in.tail11, beta=beta.assumed, rare.maf=1e-4, y.prime=y.prime, prs=prs.prime ))


#fit <- summary(lm( yyy ~ I(prs1) ))
#r <- sqrt(fit$r.squared)
#ptr <- which( rare.effects1==0 & rare.effects2==0 )
#ptr <- 1:n
#y.prime <- yyy[ptr]
#prs.prime <- prs1[ptr]

#fit <- summary( lm(prs.prime ~ y.prime ))
#print(fit)
#print( sqrt(fit$r.squared) / sd(y.prime) )
#print( sqrt(fit$r.squared) * mean(y.prime) / sd(y.prime) )

#K = 1
#mu.y = mean(y.prime)
#sigma.y = sd(y.prime)
#K.prime = (K - mu.y)/sigma.y
#r <- sqrt(fit$r.squared)
#print( mean( prs.prime[ y.prime>K ] ) )
#print( r * mean( y.prime[ y.prime>K ] ) / sigma.y - r*mu.y / sigma.y )
#print( r * dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE) )

#mu.y = 0
#K.prime = (K - mu.y)/sigma.y
#print( r * dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE) )

#mu.y = mean(y.prime)
#sigma.y = 1
#K.prime = (K - mu.y)/sigma.y
#print( r * dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE) )

#mu.y = 0
#sigma.y = 1
#K.prime = (K - mu.y)/sigma.y
#print( r * dnorm(K.prime)/pnorm(K.prime,lower.tail=FALSE) )

#print( r * ( mean( y.prime[y.prime>2] ) - mean(y.prime) ) )
#print( mean( prs.prime[y.prime>2] ) )

#test <- prs.test( prs1, yy )
#prs.test( prs2, yy )
#prs.test( (prs2+rare.effects1), yy )
