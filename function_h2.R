# mean trait value in tail / POPout effect size needs to be consistent with effect sizes of rare variants pushing individuals in to the tail
trait.hist <- function( m0, m1, maf, beta, r2=NULL, tail=0.01, n=100000 ){
    prs <- sqrt(r2) * rnorm(n)
    y.env <- sqrt(1-r2) * rnorm(n)
    y.low <- -beta * rbinom( n, round(2*m0), maf )
    y.high <- beta * rbinom( n, round(2*m1), maf )
#    y.low <- -beta * apply( replicate( m0, rbinom( n, 2, maf ) ), 1, sum )
#    y.high <- beta * apply( replicate( m1, rbinom( n, 2, maf ) ), 1, sum )
    yy <- prs + y.env
    y <- yy + y.low + y.high
    y <- y - mean(y)
    hist(y)
    y.norm <- quantile.normalise(y)
#    print(shapiro.test(y[1:5000]))
#    print( table(y.high-y.low) )
    print("QN, regression on all")
    print(prs.test( prs, y.norm ))# centre.reg=0.1 ))

    print("Raw, regression on all")
    print(prs.test( prs, y ))#, centre.reg=0.1 ))

    print("QN, regression on middle 80%")
    print(prs.test( prs, y.norm, centre.reg=0.1  ))# centre.reg=0.1 ))

    print("Raw, regression on middle 80%")
    print(prs.test( prs, y, centre.reg=0.1  ))#, centre.reg=0.1 ))

    print( var(y.low+y.high) / (1+var(y.low+y.high) ) )
#    print( mean( y<qnorm(tail) & y.low==-beta ) / mean( y<qnorm(tail) ) )
#    print( mean( y>qnorm(1-tail) & y.high==beta ) / mean( y>qnorm(1-tail) ) )

    ptr0 <- which( yy<qnorm(tail) )
    ptr <- which( y<qnorm(tail) )
#    print( mean(prs[ptr]) - mean(prs[ptr0]) )
    ptr0 <- which( yy>qnorm(tail,lower.tail=FALSE) )
    ptr <- which( y>qnorm(tail,lower.tail=FALSE) )
#    prs1 <- prs/sd(prs)
#    print( mean(prs1[ptr]) - mean(prs1[ptr0]) )
}

sim.check <- function( m0, m1, maf, beta, r2=NULL, tail=0.01, n=500000, verbose=FALSE ){
    prs <- sqrt(r2) * rnorm(n)
    y.env <- sqrt(1-r2) * rnorm(n)
    y.low <- -beta * rbinom( n, round(2*m0), maf )
    y.high <- beta * rbinom( n, round(2*m1), maf )
    yy <- prs + y.env
    y <- yy + y.low + y.high
    if( verbose ){
        print(table(y.low))
        print(table(y.high))
        print(skew(y))
        print(kurtosis(y))
    }
    y.norm <- quantile.normalise(y)
    test <- prs.test( prs, y.norm )
    h2 <- var(y.low+y.high) / (1+var(y.low+y.high))
    return( c(h2, as.numeric(test[,1])) )
}

var.explained <- function(data,ptr){
    ptr.PRS <- grep('PRS',colnames(data))
    fit0 <- summary(lm( data$y ~ 1, subset=ptr ))
    fit1 <- summary(lm( data$y ~ as.matrix(data[,ptr.PRS]), subset=ptr ))
    R2 <- 1 - (1-fit1$adj.r.squared) / (1-fit0$adj.r.squared)
    return(R2)
}

est.prop.in.tail <- function( effect.size, beta, r2, tail=0.01 ){
    percentile.pushed.to.tail <- vector()
    prop.in.tail <- c(0,0)
    kappa <- qnorm( tail, lower.tail=FALSE )
    for( i in 1:2 ){
        percentile.pushed.to.tail[i] <- pnorm( beta - kappa  )
        Z1 <- qnorm( percentile.pushed.to.tail[i], lower.tail=FALSE )
        Z0 <- qnorm( tail, lower.tail=FALSE )
#        m.Y99.rare <- r * (dnorm(Z1)-dnorm(Z0)) / ( percentile.pushed.to.tail - tail )
        m.Y99.rare <- sqrt(r2) * dnorm(Z1) / percentile.pushed.to.tail[i]
        m.Y99 <- sqrt(r2) * dnorm(Z0) / tail
        if( effect.size[i]>0 ){
            prop.in.tail[i] <- effect.size[i] / ( m.Y99 - m.Y99.rare )
        }
    }
    return( prop.in.tail )
}

h2.rare.big <- function( prop.in.tail,  tail=0.01, rare.maf=1e-4, beta ){
# beta -- presumed effect size of rare large effect variants
    kappa <- qnorm( tail, lower.tail=FALSE )
    percentile.pushed.to.tail <- pnorm(beta-kappa)

    sum.rare.freq <- vector()
    var.rare <- vector()
    m <- vector()
    for( i in 1:2 ){
# prop.in.tail[i]*tail = population proportion pushed to tail by rare large effect variants
        sum.rare.freq[i] <- prop.in.tail[i]*tail / (prop.in.tail[i]*tail + (1-prop.in.tail[i])*percentile.pushed.to.tail )
        m[i] <- 0.5*sum.rare.freq[i] / rare.maf # number of rare large effect variants
        var.rare[i] <- m[i] * 2*rare.maf*(1-rare.maf) * beta^2
    }

    h2 <- sum(var.rare) / (1+sum(var.rare))
    ret <- c( h2, sum.rare.freq, m, percentile.pushed.to.tail )
    names(ret) <- c( 'h2', 'sum.rare.freq1', 'sum.rare.freq2',
                    'm1', 'm2', 'percentile.pushed.to.tail' )
    return(as.list(ret))
}


est.prop.in.tail2 <- function( effect.size, beta, r2, tail=0.01, mu.y=0, sigma.y=1 ){
    r = sqrt(r2)
    percentile.pushed.to.tail <- vector()
    prop.in.tail <- c(0,0)
    K <- qnorm( tail, lower.tail=FALSE )
    K.prime = vector()
    K.prime.rare = vector()
    K.prime[1] = (-K - mu.y)/sigma.y
    K.prime[2] = (K - mu.y)/sigma.y
    K.prime.rare[1] = (-K + beta - mu.y)/sigma.y
    K.prime.rare[2] = (K - beta - mu.y)/sigma.y
    for( i in 1:2 ){
        tail = ifelse( i==1, TRUE, FALSE )
        m.Y99.rare <- r * dnorm( K.prime.rare[i] )/ pnorm( K.prime.rare[i], lower.tail=tail )
        m.Y99 <- r * dnorm( K.prime[i] ) / pnorm( K.prime[i], lower.tail=tail )
        if( effect.size[i]>0 ){
            prop.in.tail[i] <- effect.size[i] / ( m.Y99 - m.Y99.rare )
        }
    }
    return( prop.in.tail )
}

h2.rare.big2 <- function( prop.in.tail, tail=0.01, rare.maf=1e-4, beta, mu.y=0, sigma.y=1 ){
# beta -- presumed effect size of rare large effect variants
    K <- qnorm( tail, lower.tail=FALSE )
    K.prime = vector()
    K.prime.rare = vector()
    K.prime[1] = (-K - mu.y)/sigma.y
    K.prime[2] = (K - mu.y)/sigma.y
    K.prime.rare[1] = (-K + beta - mu.y)/sigma.y
    K.prime.rare[2] = (K - beta - mu.y)/sigma.y

    sum.rare.freq <- vector()
    var.rare <- vector()
    m <- vector()
    for( i in 1:2 ){
        tail = ifelse( i==1, TRUE, FALSE )
        sum.rare.freq[i] <- prop.in.tail[i]*pnorm( K.prime[i], lower.tail=tail ) /
            ( prop.in.tail[i]*pnorm( K.prime[i], lower.tail=tail ) +
              (1-prop.in.tail[i])*pnorm( K.prime.rare[i], lower.tail=tail ) )
        m[i] <- 0.5*sum.rare.freq[i] / rare.maf # number of rare large effect variants
        var.rare[i] <- m[i] * 2*rare.maf*(1-rare.maf) * beta^2
    }

    h2 <- sum(var.rare)
    ret <- c( h2, sum.rare.freq, m )
    names(ret) <- c( 'h2', 'sum.rare.freq1', 'sum.rare.freq2',
                    'm1', 'm2' )
    return(as.list(ret))
}

h2.est.emp <- function( n, effect.size, prs.r2, h2.common, beta, sd.beta=0, rare.maf=1e-4, tail=0.01, n.samples=50, iter.max=5 ){
    good.sim <- FALSE
    ret <- matrix(ncol=6,nrow=0)
    colnames(ret) <- c( 'h2', 'sum.rare.freq1', 'sum.rare.freq2', 'm1', 'm2', 'iter' )
    p.in.tail <- est.prop.in.tail2( effect.size, beta=beta, r2=prs.r2 )
    p.in.tail <- ifelse( p.in.tail<0, 0, p.in.tail )
    p.in.tail <- ifelse( p.in.tail>1, 1, p.in.tail )
    ex <- h2.rare.big2( p.in.tail, beta=beta, rare.maf=rare.maf, tail=tail )
    m1 <- ifelse( ex$m1<0, 0, round(ex$m1) )
    m2 <- ifelse( ex$m2<0, 0, round(ex$m2) )
    while( h2.common + 2 * rare.maf * (1-rare.maf) * (m1+m2) * beta^2 > 0.95 ){
        m1 <- m1/2
        m2 <- m2/2
    }
    m1.i <- m1
    m2.i <- m2
    for( i in 1:n.samples ){
#        print(i)
        m1.new <- 1e6
        m2.new <- 1e6
        iter.sum <- 0
        while( (h2.common + 2 * rare.maf * (1-rare.maf) * (m1.new+m2.new) * beta^2 > 0.95) & iter.sum < iter.max )
        {
            iter.sum <- iter.sum + 1
            sim.data <- sim.pheno( n=n, m1=m1, m2=m2, rare.maf=rare.maf, beta=beta,
                                  h2.common=h2.common, prs.r2=prs.r2, effect.size, sd.beta=sd.beta, tail=tail )
#            print(sim.data[[3]])
            mu.y = mean(sim.data[[2]][,'y.prime'])
            sigma.y = sd(sim.data[[2]][,'y.prime'])
            p.in.tail <- est.prop.in.tail.emp( effect.size, beta, r2=prs.r2,
                                              y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
            p.in.tail2 <- est.prop.in.tail2( effect.size, beta, r2=prs.r2, mu.y=mu.y, sigma.y=sigma.y )
            p.in.tail <- ifelse( ( is.nan(p.in.tail) | p.in.tail<0 | p.in.tail>1 ) , p.in.tail2, p.in.tail )
            p.in.tail <- ifelse( p.in.tail<0, 0, p.in.tail )
            p.in.tail <- ifelse( p.in.tail>1, 1, p.in.tail )
            ex <- h2.rare.big.emp( p.in.tail, beta=beta, rare.maf=rare.maf,
                                  y.prime=sim.data[[2]]$y.prime, prs=sim.data[[2]]$prs.prime )
            m1.new <- ifelse( ex$m1<0, 0, round(ex$m1) )
            m2.new <- ifelse( ex$m2<0, 0, round(ex$m2) )
            m1.new <- ifelse( is.na(ex$m1), 1e6, m1.new )
            m2.new <- ifelse( is.na(ex$m2), 1e6, m2.new )
#            print(unlist(ex))
#            print( h2.common + 2 * rare.maf * (1-rare.maf) * (m1.new+m2.new) * beta.assumed^2 )
        }
        if( iter.sum==iter.max ){
            m1 <- m1.i
            m2 <- m2.i
        }else{
            m1 <- m1.new
            m2 <- m2.new
            m1.i <- m1
            m2.i <- m2
        }
        if( (h2.common + 2 * rare.maf * (1-rare.maf) * (m1.new+m2.new) * beta^2 < 0.95) ){
            ret <- rbind( ret, c( unlist(ex), iter.sum ) )
        }
    }
    if( nrow(ret)>1 ){
        ret1 <- c( apply( ret, 2, mean ), nrow(ret) )
    }else{
        ret1 <- c( ret, nrow(ret) )
    }
    return( ret )
}

est.prop.in.tail.emp <- function( effect.size, beta, r2, y.prime, prs, tail=0.01 ){
    r <- sqrt(r2)
    kappa <- qnorm( tail, lower.tail=FALSE )
    prop.in.tail <- c(0,0)

    m.Y99 <- mean( prs[ y.prime < -kappa ] )
    m.Y99.rare <- mean( prs[ y.prime < (-kappa+beta) ] )
    m.Y99.null <- -r * dnorm( kappa ) / pnorm( kappa, lower.tail=FALSE )
    prop.in.tail[1] <- (-effect.size[1] - m.Y99.null + m.Y99 ) / ( m.Y99 - m.Y99.rare ) # negative since forced +ve effect for tail

    m.Y99 <- mean( prs[ y.prime > kappa ] )
    m.Y99.rare <- mean( prs[ y.prime > (kappa-beta) ] )
    m.Y99.null <- r * dnorm( kappa ) / pnorm( kappa, lower.tail=FALSE )
    prop.in.tail[2] <- ( effect.size[2] - m.Y99.null + m.Y99 ) / ( m.Y99 - m.Y99.rare )

    return( prop.in.tail )
}

h2.rare.big.emp <- function( prop.in.tail,  tail=0.01, rare.maf=1e-4, beta, y.prime, prs ){
# beta -- presumed effect size of rare large effect variants
    mu.y = mean(y.prime)
    sigma.y = sd(y.prime)
    K <- qnorm( tail, lower.tail=FALSE )
    K.prime = vector()
    K.prime.rare = vector()
    K.prime[1] = (-K - mu.y)/sigma.y
    K.prime[2] = (K - mu.y)/sigma.y
    K.prime.rare[1] = (-K + beta - mu.y)/sigma.y
    K.prime.rare[2] = (K - beta - mu.y)/sigma.y

    sum.rare.freq <- vector()
    var.rare <- vector()
    m <- c(0,0)

    sum.rare.freq[1] <- prop.in.tail[1]*mean( y.prime < -K ) /
        ( prop.in.tail[1]*mean( y.prime < -K ) + (1-prop.in.tail[1])*mean( y.prime<(-K+beta) ) )
    sum.rare.freq[2] <- prop.in.tail[2]*mean( y.prime > K ) /
        ( prop.in.tail[2]*mean( y.prime > K ) + (1-prop.in.tail[2])*mean( y.prime>(K-beta) ) )

#    print( c( prop.in.tail[2], mean( y.prime > K ), sum( y.prime > K )  ) )
#    print(sum.rare.freq)

    if( sum( y.prime < -K ) < 2 ){
        print('here1')
        print( sum( y.prime < -K ) )
        sum.rare.freq[1] <- prop.in.tail[1]*pnorm( K.prime[1], lower.tail=TRUE ) /
            ( prop.in.tail[1]*pnorm( K.prime[1], lower.tail=TRUE ) +
              (1-prop.in.tail[1])*pnorm( K.prime.rare[1], lower.tail=TRUE ) )
        }
    if( sum( y.prime > K ) < 2 ){
        print('here2')
        print(sum( y.prime > K ))
        sum.rare.freq[2] <- prop.in.tail[2]*pnorm( K.prime[2], lower.tail=FALSE ) /
            ( prop.in.tail[2]*pnorm( K.prime[2], lower.tail=FALSE ) +
              (1-prop.in.tail[2])*pnorm( K.prime.rare[2], lower.tail=FALSE ) )
    }
#    print(sum.rare.freq)

    for( i in 1:2 ){
        m[i] <- 0.5*sum.rare.freq[i] / rare.maf # number of rare large effect variants
        var.rare[i] <- m[i] * 2*rare.maf*(1-rare.maf) * beta^2
    }

    h2 <- sum(var.rare)
    ret <- c( h2, sum.rare.freq, m )
    names(ret) <- c( 'h2', 'sum.rare.freq1', 'sum.rare.freq2',
                    'm1', 'm2' )
    return(as.list(ret))
}

sim.pheno <- function( n, m1, m2, rare.maf, beta, h2.common, prs.r2, effect.size, sd.beta=0, tail=0.01 ){

    #    while( s.upper<5 | s.lower<5 ){
#    K <- qnorm( tail, lower.tail=FALSE )
#    s.lower <- 0
#    s.upper <- 0
#        s.lower <- sum( y.prime < -K )
#        s.upper <- sum( y.prime > K )
#        print( c( s.upper, s.lower ) )

    test <- matrix(nrow=2,ncol=5,data=0)
    colnames(test) <- c( 'effect', 'ci.lower', 'ci.upper', 'p', 'se' )
    rownames(test) <- c( 'lower', 'upper' )
    common.effects <- sqrt(h2.common) * rnorm( n=n )
    if( m1>0 ){
        rare.effects1 <- apply( matrix(ncol=n,data=replicate( n, rnorm( m1, -beta, sd.beta ) *
                                                                 rbinom(n=m1,p=rare.maf,size=2))), 2, sum )
    }else{
        rare.effects1 <- rep(0,n)
    }
    if( m2>0 ){
        rare.effects2 <- apply( matrix(ncol=n,data=replicate( n, rnorm( m2, beta, sd.beta ) *
                                                                 rbinom(n=m2,p=rare.maf,size=2))), 2, sum )
    }else{
        rare.effects2 <- rep(0,n)
    }
    h2.env <- 1 - h2.common - var(rare.effects1) - var(rare.effects2)
    env <- sqrt(h2.env) * rnorm( n=n )
    y <- common.effects + rare.effects1 + rare.effects2 + env
    yy <- quantile.normalise(y)
    prs.e <- h2.common * ( h2.common / prs.r2 - 1 )
    prs <- common.effects + sqrt(prs.e) * rnorm( n=n )
    delta.prs.r2 <- cor( prs, yy )^2 - prs.r2
    prs.ee <- prs.e / 10
    while( delta.prs.r2>0.001 ){
        prs <- prs + sqrt(prs.ee) * rnorm( n=n )
        delta.prs.r2 <- cor( prs, yy )^2 - prs.r2
    }
    prs <- prs / sd(prs)
    test <- prs.test( prs, yy )
#    print(test)
    good.sim <-
        ifelse( test['lower','ci.lower'] < effect.size[1]&effect.size[1] < test['lower','ci.upper'] &
                test['upper','ci.lower'] < effect.size[2]&effect.size[2] < test['upper','ci.upper'],
               TRUE, FALSE )

    ptr <- which( rare.effects1==0 & rare.effects2==0 )
    y.prime <- yy[ptr]
    prs.prime <- prs[ptr]

    return( list( data.frame( y, yy, prs, rare.effects1, rare.effects2 ),
                 data.frame( y.prime, prs.prime ),good.sim ) )
}

sim <- function( n, m, rare.maf, beta, prs.r2, centre.reg=NULL ){
    y <- sim.pheno( n, m, rare.maf, beta, prs.r2 )[,1]
    h2 <- var(sum.rare) / var(y)

    qq <- quantile( y, c( 0.01, 0.99 ) )
    print( length(which( y<qq[1] & rare.effects[,1] )>0 ) / sum(y<qq[1]) )
    print( length(which( y>qq[2] & rare.effects[,2] )>0 ) / sum(y>qq[2]) )
    prs.test.plot( prs1=prs/sd(prs), y=y )

    test <- prs.test( prs1=prs/sd(prs), y=y )#, centre.reg=centre.reg )
    ret <- list( h2, test )
    names(ret) <- c( 'h2', 'test' )

    return( ret )
}

est.prop.in.tail.old <- function( effect.size, beta, r, tail=0.01 ){
    kappa <- qnorm( tail, lower.tail=FALSE )
    percentile.pushed.to.tail <- pnorm(beta-kappa)
    prop.in.tail <- c(0,0)
    for( i in 1:2 ){
        a.rare <- qnorm( percentile.pushed.to.tail, lower.tail=FALSE )
        m.rare <- r * dnorm(a.rare) / percentile.pushed.to.tail
        a.prs <- qnorm(tail,lower.tail=FALSE)
        m.prs <- r * dnorm(a.prs) / tail
        if( effect.size[i]>0 ){
            prop.in.tail[i] <- effect.size[i] / ( m.prs - m.rare )
        }
    }
    return( prop.in.tail )
}

h2.rare.big.old <- function( prop.in.tail,  tail=0.01, rare.maf=1e-4, beta ){
# beta -- presumed effect size of rare large effect variants
    kappa <- qnorm( tail, lower.tail=FALSE )
    percentile.pushed.to.tail <- pnorm(beta-kappa)

    sum.rare.freq <- vector()
    var.rare <- vector()
    m <- vector()
    for( i in 1:2 ){
# prop.in.tail[i]*tail = population proportion pushed to tail by rare large effect variants
        sum.rare.freq[i] <- prop.in.tail[i]*tail / percentile.pushed.to.tail
        m[i] <- 0.5*sum.rare.freq[i] / rare.maf # number of rare large effect variants
        var.rare[i] <- m[i] * 2*rare.maf*(1-rare.maf) * beta^2
    }

    h2 <- sum(var.rare) / (1+sum(var.rare))
    ret <- c( h2, sum.rare.freq, m, percentile.pushed.to.tail )
    names(ret) <- c( 'h2', 'sum.rare.freq1', 'sum.rare.freq2',
                    'm1', 'm2', 'percentile.pushed.to.tail' )
    return(as.list(ret))
}
