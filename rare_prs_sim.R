source('~/bin/functions.R')
source('functions.R')
source('function_h2.R')
n <- 2e5
h2.common <- 0.3
prs.r2.1 <- 0.15
prs.r2.2 <- 0.1
prs1.e <- h2.common * ( h2.common / prs.r2.1 - 1 )
prs2.e <- h2.common * ( h2.common / prs.r2.2 - 1 )
beta.real <- 2
beta.assumed <- 2

symm <- TRUE
iter <- 1
h2.est1 <- matrix( nrow=iter, ncol=5 )

rare.maf <- 1e-4
n.rare.effects1 <- 40
n.rare.effects2 <- 30

for( i in 1:iter ){
    common.effects <- sqrt(h2.common) * rnorm( n=n )
    rare.effects1 <- apply(matrix(ncol=n,data=replicate( n,
                                                        rnorm( n.rare.effects1, -beta.real, 1 ) *
                                                        rbinom(n=n.rare.effects1,p=rare.maf,size=2))), 2, sum )
    rare.effects2 <- apply(matrix(ncol=n,data=replicate( n,
                                                        rnorm( n.rare.effects2, beta.real, 1 ) *
                                                        rbinom(n=n.rare.effects2,p=rare.maf,size=2))), 2, sum )
#    if( symm ){
#        rare.effects1 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects1 # captured rare effects
#        rare.effects2 <- (2*rbinom( n=n, p=0.99, size=1 ) - 1) * rare.effects2 # missed rare effects
#    }

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
    print(test)
    p.in.tail1 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex1 <- h2.rare.big2( p.in.tail1, beta=beta.assumed, rare.maf=1e-4, mu.y=mu.y, sigma.y=sigma.y )

    fit <- summary(lm( yyy ~ I(prs2) ))
    test <- prs.test( (prs2), yyy )
    p.in.tail2 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
    ex2 <- h2.rare.big2( p.in.tail2, beta=beta.assumed, mu.y=mu.y, sigma.y=sigma.y )

    fit <- summary(lm( yyy ~ I(prs1) ))
    test <- prs.test( (prs1), yyy )
    ex3 <- as.list( h2.est.emp( n=n, effect.size=test[,'effect'], prs.r2=fit$r.squared, h2.common=h2.common,
                      beta=beta.assumed, sd.beta=0 ) )

    fit <- summary(lm( yyy ~ I(prs2) ))
    test <- prs.test( (prs2), yyy )
    ex4 <- as.list( h2.est.emp( n=n, effect.size=test[,'effect'], prs.r2=fit$r.squared, h2.common=h2.common,
                      beta=beta.assumed, sd.beta=0 ) )

#    print( c( p.in.tail1, p.in.tail2, p.in.tail3, p.in.tail4 ) )
    print( c( 1-var(yy[ptr]), ex1$h2, ex2$h2, ex3$h2, ex4$h2 ) )
    h2.est1[i,] <- c( 1-var(yy[ptr]), ex1$h2, ex2$h2, ex3$h2, ex4$h2 )
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


mu.y=0
sigma.y=1

delta <- 0
fit <- summary(lm( yyy ~ I(prs1) ))
test <- prs.test( (prs1), yyy )
p.in.tail3 <- est.prop.in.tail2( test[,'effect'], (beta.assumed-delta), r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
ex3 <- h2.rare.big2( p.in.tail3, beta=(beta.assumed-delta), mu.y=mu.y, sigma.y=sigma.y  )
sim.data <- sim.pheno( n=n, m1=ex3$m1, m2=ex3$m2, rare.maf=1e-4, beta=(beta.assumed-delta),
                      h2.common=h2.common, prs.r2=fit$r.squared, sd.beta=0 )
for( i in 1:10 ){
    p.in.tail3 <- est.prop.in.tail.emp( test[,'effect'], (beta.assumed-delta),
                                       y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
    ex3 <- h2.rare.big.emp( p.in.tail3, beta=(beta.assumed-delta),
                           y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
    sim.data <- sim.pheno( n=n, m1=ex3$m1, m2=ex3$m2, rare.maf=1e-4, beta=(beta.assumed-delta),
                          h2.common=h2.common, prs.r2=fit$r.squared, sd.beta=1 )
    mu.y = mean(sim.data[[2]][,'y.prime'])
    sigma.y = sd(sim.data[[2]][,'y.prime'])
    fit.sim <- summary(lm( sim.data[[1]][,'yy'] ~ sim.data[[1]][,'prs'] ))
    test.sim <- prs.test( sim.data[[1]][,'prs'], sim.data[[1]][,'yy'] )
    print( signif( c( mu.y, sigma.y, ex3$h2, test.sim$effect, test.sim$p ), 3 ) )
}

#    mu.y=0
#    sigma.y=1
#    fit <- summary(lm( yyy ~ I(prs2) ))
#    test <- prs.test( (prs2), yyy )
#    p.in.tail4 <- est.prop.in.tail2( test[,'effect'], beta.assumed, r2=fit$r.squared, mu.y=mu.y, sigma.y=sigma.y )
#    ex4 <- h2.rare.big2( p.in.tail4, beta=beta.assumed, mu.y=mu.y, sigma.y=sigma.y  )
#    sim.data <- sim.pheno( n=n, m1=ex4$m1, m2=ex4$m2, rare.maf=1e-4, beta=beta.assumed,
#                          h2.common=h2.common, h2.env==1-h2.common, prs.r2=fit$r.squared, sd.beta=0 )
#    for( i in 1:15 ){
#        p.in.tail4 <- est.prop.in.tail.emp( test[,'effect'], beta.assumed,
#                                           y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
#        ex4 <- h2.rare.big.emp( p.in.tail4, beta=beta.assumed,
#                            y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
#        sim.data <- sim.pheno( n=n, m1=ex4$m1, m2=ex4$m2, rare.maf=1e-4, beta=beta.assumed,
#                              h2.common=h2.common, h2.env==1-h2.common, prs.r2=fit$r.squared, sd.beta=0 )
#        mu.y = mean(sim.data[[2]][,'y.prime'])
#        sigma.y = sd(sim.data[[2]][,'y.prime'])
#        fit.sim <- summary(lm( sim.data[[1]][,'yy'] ~ sim.data[[1]][,'prs'] ))
#        test.sim <- prs.test( sim.data[[1]][,'prs'], sim.data[[1]][,'yy'] )
#        print( signif( c( mu.y, sigma.y, ex4$h2, test.sim$effect, test.sim$p ), 3 ) )
#    }
