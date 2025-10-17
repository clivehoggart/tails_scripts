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

est.prop.in.tail <- function( effect.size, beta, r2, tail=0.01, mu=0 ){
    percentile.pushed.to.tail <- vector()
    prop.in.tail <- c(0,0)
    kappa <- qnorm( tail, lower.tail=FALSE )
    sign <- 1
    for( i in 1:2 ){
        sign <- -sign
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

est.prop.in.tail <- function( effect.size, beta, r2, tail=0.01, mu.y=0, sigma.y=1 ){
    percentile.pushed.to.tail <- vector()
    prop.in.tail <- c(0,0)
    K <- qnorm( tail, lower.tail=FALSE )
    K.prime[1] = (-K - mu.y)/sigma.y
    K.prime[2] = (K - mu.y)/sigma.y
    K.prime.rare[1] = (-K + beta - mu.y)/sigma.y
    K.prime.rare[2] = (K - beta - mu.y)/sigma.y
    tail = FALSE
    for( i in 1:2 ){
        m.Y99.rare <- sqrt(r2) * (mu.y + sigma.y*dnorm(K.prime.rare)/pnorm(K.prime.rare,lower.tail=tail)) / sigma.y - mu.y / sigma.y )
        m.Y99 <- sqrt(r2) * (mu.y + sigma.y*dnorm(K.prime)/pnorm(K.prime,lower.tail=tail)) / sigma.y - mu.y / sigma.y )
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

sim.pheno <- function( n, m, rare.maf, beta, prs.r2 ){
    rare.genos <- matrix( ncol=2, nrow=n, data=0 )
    rare.effects <- matrix( ncol=2, nrow=n, data=0 )
    for( i in 1:2 ){
        if( m[i]>0 ){
            rare.genos[,i] <- apply( replicate( m[i], rbinom( n, 2, rare.maf ) ), 1, sum )
            rare.effects[,i] <- rare.genos[,i] * beta * (-1)^i
        }
    }
    sum.rare <- apply( rare.effects, 1, sum )
    prs <- rnorm(n)*sqrt(prs.r2)
    y <- rnorm(n)*sqrt(1-prs.r2-var(sum.rare)) + prs + sum.rare
    y <- y - mean(y)
    return( cbind( y, sum.rare ) )
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

est.prop.in.tail2 <- function( effect.size, beta, r, tail=0.01 ){
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

h2.rare.big2 <- function( prop.in.tail,  tail=0.01, rare.maf=1e-4, beta ){
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
