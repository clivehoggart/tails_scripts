library(MASS)
library(parallel)
library(data.table)

library(lmtest)

assign.to.quants <- function( y, n ){
    q <- quantile( y, 1:(n-1)/n, na.rm=TRUE )
    y.cat <- sapply( lapply( y, '>=', q ), sum )
    return(y.cat)
}

loglike.prs.r2 <- function( r2, prs, y ){
    n <- length(y)
    sigma2 <- r2*(1-r2)
    l <- n*log(sigma2) + sum( (prs - y*r2)^2 ) / sigma2
    return(l)
}

prs.test1 <- function( prs, y, K, K1, tail='upper' ){#heteroskedasticity test
    if( tail=='upper' ){
        T <- quantile( y, 1-K )
        ptr <- y>T
    }
    if( tail=='lower' ){
        T <- quantile( y, K )
        ptr <- y<T
    }
    if( tail=='quantile' ){
        T1 <- quantile( y, K )
        T2 <- quantile( y, K1 )
        ptr <- T1 < y&y < T2
    }
    print(sum(ptr))

    p <- bptest( prs~y, varformula=~I(ptr) )$p.value

    return(p)
}

prs.test <- function( prs1, y, K=0.01, tail="two" ){# Test for residuals mean=0 in the tail(s)
    prs1 <- prs1/sd(prs1)
    fit <- lm( prs1 ~ y )

    # Test upper tail "regression to mean"
    T <- quantile( y, 1-K, na.rm=TRUE )
    ptr <- which(y>T)
    alt <- ifelse( tail=="one", "less", tail )
    test <- t.test( fit$residual[ptr], alternative=alt )
    p2 <- test$p.value
    effect2 <- -test$estimate
    ci2 <- c( -test$conf.int[2], -test$conf.int[1] )
    se2 <- test$stderr

    # Test lower tail "regression to mean"
    T <- quantile( y, K, na.rm=TRUE )
    ptr <- which(y<T)
    alt <- ifelse( tail=="one", "greater", tail )
    test <- t.test( fit$residual[ptr], alternative=alt )
    p1 <- test$p.value
    effect1 <- test$estimate
    ci1 <- test$conf.int
    se1 <- test$stderr

    ret <- as.data.frame(rbind( c( effect1, ci1, p1, se1), c( effect2, ci2, p2, se2 ) ) )
    colnames(ret) <- c( 'effect', 'ci.lower', 'ci.upper', 'p', 'se' )
    rownames(ret) <- c( 'lower', 'upper' )

    return(ret)
}

prs.test.plot <- function( prs1, y ){
    prs1 <- prs1/sd(prs1)
    fit <- lm( prs1 ~ y )

    q <- assign.to.quants( y, 100 )
    m.o <- tapply( prs1, q, mean )
    ptr.use <- setdiff( 1:length(y), fit$na.action )
    m.e <- tapply( fit$fitted.values, q[ptr.use], mean )
    plot( 1:100, m.o, xlab='Trait quantile', ylab='Mean PRS', ylim=range(m.e) )
    lines( 1:100, m.e, col='red' )
}

prs.test.qc <- function( prs, y, n.q=100, K=0.1, K1=0.9 ){
    q <- assign.to.quants( y, n.q )
    ptr <- which( n.q*K <= q&q < n.q*K1 )
    fit <- lm( prs ~ y, subset=ptr )
    p <- vector()
    for( j in sort(unique(q[ptr])) ){
        ptr1 <- which( q[ptr]==j )
#        print(range(y[ptr[ptr1]]))
        p <- c( p, t.test( fit$residuals[ptr1] )$p.value )
    }
    qc.test <- ifelse( min(p) < 0.05/length(unique(q[ptr])), 'FAIL', 'PASS' )
    return(c(qc.test,min(p)))
}

prs.test.emp <- function( prs, y, n.q=100 ){
    stat <- vector()
    fit <- lm( prs ~ y )

    q <- assign.to.quants( y, n.q )
    for( i in 1:n.q ){
        ptr <- which(q==(i-1))
        test <- t.test( fit$residual[ptr] )
        stat[i] <- test$statistic
    }
    r1 <- rank(-stat)
    r2 <- rank(stat)
    p <- c( r1[1]/n.q, r2[n.q]/n.q )
    return(p)
}

prs.test3 <- function( prs, y, K, tail='upper' ){# Motivated by sibs test, unnescarrily complicated
    r2 <- optim( 0.5, loglike.prs.r2, prs=prs, y=y,
                method='Brent', lower=0, upper=1  )$par
    T <- qnorm( K, lower.tail=tail=='lower' )

# Moments of truncated trait distribution -- y
    m.trunc <- (dnorm(T) / K) * (-1)^(tail=='lower')
    v.trunc  <- (1 + T*m.trunc - m.trunc*m.trunc)

    if( tail=='upper' )
        ptr <- which(y>T)
    if( tail=='lower' )
        ptr <- which(y<T)

# Moments of truncated PRS distribution
    m.trunc.prs <- m.trunc*r2
    v.trunc.prs <- v.trunc*r2*r2 + r2*(1-r2)

# Score test
    z <- (prs[ptr] - m.trunc.prs) / v.trunc.prs
    nn <- length(ptr)
    U <- sum(z)
    I <- -nn / v.trunc.prs
    z.stat <- U / sqrt(-I)
    p <- pnorm( z.stat, lower.tail=tail=='upper' )
    return(p)
}

loglike <- function( h2, s ){
    sigma2 <- 1-h2^2/4
    n <- nrow(s)
    f <- sum( (s[,2] - s[,1]*h2/2)^2 ) / (sigma2) + n*log(sigma2)
    return(f)
}

h2.est <- function( h2.init=0.5, s ){
    fit <- optim( h2.init, loglike, s=s, lower=0, upper=1, method="L-BFGS-B" )
    return(fit$par)
}

denovo.test <- function( s, h2, K ){
    T <- c( qnorm(K), -qnorm(K) ) # Gaussian truncation point
    z <- (s[,2] - s[,1]*h2/2) / (1-h2^2/4)
    ptr <- list( which( s[,1]<T[1] ), which( s[,1]>T[2] ) )
    p <- vector()
    for( i in 1:2 ){
        n <- length(ptr[[i]])
        U <- sum(z[ptr[[i]]])
        I <- -n / (1-h2^2/4)
        z.stat <- U / sqrt(-I)
        p[i] <- pnorm( z.stat, lower.tail=FALSE )
        p[i] <- ifelse( i==1, p[i], 1-p[i] )
    }
    return(p)
}

mend.test <- function( s, h2, K ){
    T <- -qnorm(K) # =2 Gaussian truncation point
    i <- dnorm(T)/K # mean of normal truncated at T
    z.affected <- (T-i*h2/2) / sqrt( 1-h2^2*i*(i-T)/4 )
    p.affected <- pnorm( z.affected, lower.tail=FALSE )

    p.con <- vector()
    for( i in 1:2 ){
        ii <- (-1)^i
        ptr1 <- which( ii * s[,1] > T )
        ptr2 <- which( ii * s[,2] > T )
        observed.concordance <- length(intersect(ptr1,ptr2))
        p.con[i] <- binom.test( observed.concordance, length(ptr1), p.affected,
                            alternative="greater" )$p.value
    }

    return(p.con)
}

