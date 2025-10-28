library(data.table)
source('~/bin/functions.R')
source('~/tails_scripts/functions.R')
source('~/tails_scripts/function_h2.R')
x <- fread('~/tails/results/tade_exome_results.txt')
x.common <- fread('~/tails/results/CLIVE_REQUEST.txt',sep=' ')
x.se <- fread('~/tails/results/SE_MATCH.txt')
colnames(x.se)[c(1,4,5,6,7)] <- c("trait","lower.se.common", "upper.se.common","lower.se.combo", "upper.se.combo")
traits <- unique(x$V1)
iter <- 1000
beta <- c( 2, 3, 4 )
rare.h2 <- array( dim=c( length(traits), length(beta), iter ) )
rare.lower <- array( dim=c( length(traits), length(beta), iter ) )
rare.upper <- array( dim=c( length(traits), length(beta), iter ) )
theta.lower <- array( dim=c( length(traits), length(beta), iter ) )
theta.upper <- array( dim=c( length(traits), length(beta), iter ) )
p.in.tail.upper <- array( dim=c( length(traits), length(beta), iter ) )
p.in.tail.lower <- array( dim=c( length(traits), length(beta), iter ) )

x.common <- x.common[match(traits,x.common$V1),]
just.common <- TRUE

for( i in 1:length(traits) ){
    if( length(which( x$V1==traits[i] ))==9 ){
        tmp <- x$V5[which( x$V1==traits[i] & x$V2=="comboResult" )]
        r2 <- as.numeric(strsplit(tmp,',')[[1]][1])
        h2 <- as.numeric(strsplit( x.common$V3[i], '=' )[[1]][2])

        ptr.se <- match( traits[i], x.se$trait )
        if( !just.common ){
            tmp <- x$V4[which( x$V1==traits[i] & x$V2=="comboTail:1%" )]
            lower.effect <- as.numeric(strsplit(tmp,',')[[1]][1])
            tmp <- x$V4[which( x$V1==traits[i] & x$V2=="comboTail:99%" )]
            upper.effect <- as.numeric(strsplit(tmp,',')[[1]][1])
            lower.effect <- ifelse( lower.effect<0, 0, lower.effect )
            upper.effect <- ifelse( upper.effect<0, 0, upper.effect )
            lower.se <- x.se$lower.se.combo[ptr.se]
            upper.se <- x.se$upper.se.combo[ptr.se]
        }else{
            lower.effect <- as.numeric(strsplit( x.common$V4[i], ',' )[[1]][2])
            upper.effect <- as.numeric(strsplit( x.common$V5[i], ',' )[[1]][2])
            lower.se <- x.se$lower.se.common[ptr.se]
            upper.se <- x.se$upper.se.common[ptr.se]
        }

        combo.stats <- c( r2, lower.effect, upper.effect, lower.se, upper.se )
        lower.effect.sample <- rnorm( n=iter, mean=lower.effect, sd=lower.se )
        upper.effect.sample <- rnorm( n=iter, mean=upper.effect, sd=upper.se )

        lower.effect.sample <- ifelse( lower.effect.sample<0, 0, lower.effect.sample )
        upper.effect.sample <- ifelse( upper.effect.sample<0, 0, upper.effect.sample )

        for( k in 1:iter ){
            for( j in 1:length(beta) ){
#                p.in.tail <- est.prop.in.tail( c( lower.effect.sample[k], upper.effect.sample[k] ), beta=beta[j], r2=r2 )
#                p.in.tail.lower[i,j,k] <- p.in.tail[1]
#                p.in.tail.upper[i,j,k] <- p.in.tail[2]
#                p.in.tail <- ifelse( p.in.tail>1, 1, p.in.tail )
#                ex <- h2.rare.big( p.in.tail, beta=beta[j], rare.maf=1e-5 )
                ex <- h2.est.emp( n=1e5, effect.size=c( lower.effect.sample[k], upper.effect.sample[k] ),
                           beta=beta[j], prs.r2=r2, h2.common=h2 )
                rare.h2[i,j,k] <- ex$h2
                rare.lower[i,j,k] <- ex$m1
                rare.upper[i,j,k] <- ex$m2
#                theta.lower[i,j,k] <- p.in.tail[1]
#                theta.upper[i,j,k] <- p.in.tail[2]
            }
        }
    }
}

out20 <- matrix(ncol=13,nrow=length(traits))
out30 <- matrix(ncol=13,nrow=length(traits))
out40 <- matrix(ncol=13,nrow=length(traits))
out10 <- matrix(ncol=13,nrow=length(traits))
out05 <- matrix(ncol=13,nrow=length(traits))
for( i in 1:length(traits) ){
    if( length(which( x$V1==traits[i] ))==9 ){
        out20[i,] <- c( round( quantile( rare.h2[i,1,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.lower[i,1,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.upper[i,1,], c(0.5,0.025,0.975) ), 3),
                       round( c( median(rare.lower[i,1,]), median(rare.upper[i,1,]) ), 1 ),
                       round( c( mean(p.in.tail.lower[i,1,]<1),
                                mean(p.in.tail.upper[i,1,]<1) ), 2 ))

        out30[i,] <- c( round( quantile( rare.h2[i,2,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.lower[i,2,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.upper[i,2,], c(0.5,0.025,0.975) ), 3),
                       round( c( median(rare.lower[i,2,]), median(rare.upper[i,2,]) ), 1 ),
                       round( c( mean(p.in.tail.lower[i,2,]<1),
                                mean(p.in.tail.upper[i,2,]<1) ), 2 ))

        out40[i,] <- c( round( quantile( rare.h2[i,3,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.lower[i,3,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.upper[i,3,], c(0.5,0.025,0.975) ), 3),
                       round( c( median(rare.lower[i,3,]), median(rare.upper[i,3,]) ), 1 ),
                       round( c( mean(p.in.tail.lower[i,3,]<1),
                                mean(p.in.tail.upper[i,3,]<1) ), 2 ))

        out10[i,] <- c( round( quantile( rare.h2[i,4,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.lower[i,4,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.upper[i,4,], c(0.5,0.025,0.975) ), 3),
                       round( c( median(rare.lower[i,4,]), median(rare.upper[i,4,]) ), 1 ),
                       round( c( mean(p.in.tail.lower[i,4,]<1),
                                mean(p.in.tail.upper[i,4,]<1) ), 2 ))

        out05[i,] <- c( round( quantile( rare.h2[i,5,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.lower[i,5,], c(0.5,0.025,0.975) ), 3),
                       round( quantile( theta.upper[i,5,], c(0.5,0.025,0.975) ), 3),
                       round( c( median(rare.lower[i,5,]), median(rare.upper[i,5,]) ), 1 ),
                       round( c( mean(p.in.tail.lower[i,5,]<1),
                                mean(p.in.tail.upper[i,5,]<1) ), 2 ))
    }
}
out05 <- data.frame( traits, out05 )
out10 <- data.frame( traits, out10 )
out20 <- data.frame( traits, out20 )
out30 <- data.frame( traits, out30 )
out40 <- data.frame( traits, out40 )
colnames(out05) <- c('trait','h2','h2.2.5%','h2.97.5%',
                     'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
                     'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
                     'm.lt','m.ut','val.lt','val.ut')
colnames(out10) <- c('trait','h2','h2.2.5%','h2.97.5%',
                     'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
                     'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
                     'm.lt','m.ut','val.lt','val.ut')
colnames(out20) <- c('trait','h2','h2.2.5%','h2.97.5%',
                     'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
                     'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
                     'm.lt','m.ut','val.lt','val.ut')
colnames(out30) <- c('trait','h2','h2.2.5%','h2.97.5%',
                     'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
                     'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
                     'm.lt','m.ut','val.lt','val.ut')
colnames(out40) <- c('trait','h2','h2.2.5%','h2.97.5%',
                     'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
                     'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
                     'm.lt','m.ut','val.lt','val.ut')

if( !just.common ){
    write.table( out05, "h2_rare_beta05.txt", quote=F, row.names=F )
    write.table( out10, "h2_rare_beta1.txt", quote=F, row.names=F )
    write.table( out20, "h2_rare_beta2.txt", quote=F, row.names=F )
    write.table( out30, "h2_rare_beta3.txt", quote=F, row.names=F )
    write.table( out40, "h2_rare_beta4.txt", quote=F, row.names=F )
}else{
    write.table( out05, "h2_rare_beta05_justcommon.txt", quote=F, row.names=F )
    write.table( out10, "h2_rare_beta1_justcommon.txt", quote=F, row.names=F )
    write.table( out20, "h2_rare_beta2_justcommon.txt", quote=F, row.names=F )
    write.table( out30, "h2_rare_beta3_justcommon.txt", quote=F, row.names=F )
    write.table( out40, "h2_rare_beta4_justcommon.txt", quote=F, row.names=F )
}
