#' test three factor interaction
#'
#'Using Parametric Bootstrap to simulate a distribution and find a p-value for the test
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@usage alg.ABC(ns, ybars, s2, a, b, c, L)
#' @param ns	sample size for each group
#' @param ybars	sample mean for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return Q:	p_value for the three factor interaction test
#'
#'@examples
#'
#'#See Q.ABmc
#'
#' @export
alg.ABC <- function(ns, ybars, s2, a, b, c, L){
  S <- diag(s2/ns) ##make S matrix

  ##make terms for X matrix
  J.abc <- rep(1, a*b*c)
  I.a <- diag(a)
  I.b <- diag(b)
  I.c <- diag(c)
  J.bc <- rep(1, b*c)
  J.a <- rep(1, a)
  J.b <- rep(1,b)
  J.c <- rep(1,c)
  I.ab <- diag(a*b)
  I.bc <- diag(b*c)

  X <- as.matrix(cbind(
    J.abc, kronecker(I.a, J.bc), kronecker(J.a, kronecker(I.b, J.c)), kronecker(J.a, kronecker(J.b, I.c)), kronecker(I.ab, J.c), kronecker(I.a, kronecker(J.b, I.c)), kronecker(J.a, I.bc)))

  #test statistic
 # library(MASS)
  SI <- t(ybars)%*%solve(S)%*%ybars -
    t(ybars)%*%solve(S)%*%X%*%ginv(t(X)%*%solve(S)%*%X)%*%t(X)%*%solve(S)%*%ybars

  ##Q, counts how many times test stat is less than PB pivot variable
  Q <- NULL
  for(j in 1:L) {
    ybar.B <- NULL
    S2B <-NULL
    for (i in 1:length(ybars)) {
      ybar.B[i] <- rnorm(1, mean=0, sd=sqrt(s2/ns)[i])  ##create bootstrap mean vector
      S2B[i] <- rchisq(1, df=(ns[i]-1)) * s2[i]/(ns[i]-1) ##create bootstrap variances vector
    }
    SB <- diag(S2B/ns)

    ##PB variable:
    SIB <-  t(ybar.B)%*%solve(SB)%*%ybar.B -
      t(ybar.B)%*%solve(SB)%*%X%*%ginv(t(X)%*%solve(SB)%*%X)%*%t(X)%*%solve(SB)%*%ybar.B

    Q[j] <- ifelse(SIB>SI, 1, 0)
  }
  return(sum(Q)/length(Q))  ##p-value
}



#' test two factor interaction
#'
#'Using Parametric Bootstrap to simulate a distribution and find a p-value for the test.
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@usage alg.BC(ns, ybars, s2, a, b, c, L)
#'
#' @param ns	sample size for each group
#' @param ybars	sample mean for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return Q:	p_value for the two factor interaction test
#'
#'@examples
#'
#'#See Q.ABmc
#'# note that the ns, ybars and s2 vectors need to be in the order reflecting subscripts
#'
#'# 111, 112, 113...., 121, 122, 123, ... , ... abc. The summarySE function from the package
#'
#'# Rmisc is handy for doing this.  The order the user specifies the "groupvars" argument will
#'
#'# put the factors in order A, B, C.  This order will matter when testing different two-way
#'
#'# interaction terms and different main effects.  See comments in the potato example.
#'
#' @export
alg.BC <- function(ns, ybars, s2, a, b, c, L){
  S <- diag(s2/ns) ##make S matrix

  ##make terms for X matrix
  J.abc <- rep(1, a*b*c)
  I.a <- diag(a)
  I.b <- diag(b)
  I.c <- diag(c)
  J.bc <- rep(1, b*c)
  J.a <- rep(1, a)
  J.b <- rep(1,b)
  J.c <- rep(1,c)
  I.ab <- diag(a*b)
  I.bc <- diag(b*c)

  X <- as.matrix(cbind(
    J.abc, kronecker(I.a, J.bc), kronecker(J.a, kronecker(I.b, J.c)), kronecker(J.a, kronecker(J.b, I.c)),
    kronecker(I.ab, J.c), kronecker(I.a, kronecker(J.b, I.c))))

  #test statistic
 # library(MASS)
  SI <- t(ybars)%*%solve(S)%*%ybars -
    t(ybars)%*%solve(S)%*%X%*%ginv(t(X)%*%solve(S)%*%X)%*%t(X)%*%solve(S)%*%ybars

  ##Q, counts how many times test stat is less than PB pivot variable
  Q <- NULL
  for(j in 1:L) {
    ybar.B <- NULL
    S2B <-NULL
    for (i in 1:length(ybars)) {
      ybar.B[i] <- rnorm(1, mean=0, sd=sqrt(s2/ns)[i])  ##create bootstrap mean vector
      S2B[i] <- rchisq(1, df=(ns[i]-1)) * s2[i]/(ns[i]-1) ##create bootstrap variances vector
    }
    SB <- diag(S2B/ns)

    ##PB variable:
    SIB <-  t(ybar.B)%*%solve(SB)%*%ybar.B -
      t(ybar.B)%*%solve(SB)%*%X%*%ginv(t(X)%*%solve(SB)%*%X)%*%t(X)%*%solve(SB)%*%ybar.B

    Q[j] <- ifelse(SIB>SI, 1, 0)
  }
  return(sum(Q)/length(Q))  ##p-value
}



#' tests factor C when only AB interaction is present.
#'
#'Using Parametric Bootstrap to simulate a distribution and find a p-value for the test
#'
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom Rmisc summarySE
#'@usage alg.C.AB(ns, ybars, s2, a, b, c, L)
#' @param ns	sample size for each group
#' @param ybars	sample mean for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return Q:	p_value for the test for factor C main effect when only AB interaction is present.
#'
#'@examples
#'
#'#See Q.ABmc
#'
#'
#'
#'
#' @export
alg.C.AB <- function(ns, ybars, s2, a, b, c, L){
  S <- diag(s2/ns) ##make S matrix

  ##make terms for X matrix
  J.abc <- rep(1, a*b*c)
  I.a <- diag(a)
  I.b <- diag(b)
  I.c <- diag(c)
  J.bc <- rep(1, b*c)
  J.a <- rep(1, a)
  J.b <- rep(1,b)
  J.c <- rep(1,c)
  I.ab <- diag(a*b)
  I.bc <- diag(b*c)

  X <- as.matrix(cbind( J.abc, kronecker(I.a, J.bc), kronecker(J.a, kronecker(I.b, J.c)), kronecker(I.ab, J.c)  ))

  #test statistic
 # library(MASS)
  SI <- t(ybars)%*%solve(S)%*%ybars -
    t(ybars)%*%solve(S)%*%X%*%ginv(t(X)%*%solve(S)%*%X)%*%t(X)%*%solve(S)%*%ybars

  ##Q, counts how many times test stat is less than PB pivot variable
  Q <- NULL
  for(j in 1:L) {
    ybar.B <- NULL
    S2B <-NULL
    for (i in 1:length(ybars)) {
      ybar.B[i] <- rnorm(1, mean=0, sd=sqrt(s2/ns)[i])  ##create bootstrap mean vector
      S2B[i] <- rchisq(1, df=(ns[i]-1)) * s2[i]/(ns[i]-1) ##create bootstrap variances vector
    }
    SB <- diag(S2B/ns)

    ##PB variable:
    SIB <-  t(ybar.B)%*%solve(SB)%*%ybar.B -
      t(ybar.B)%*%solve(SB)%*%X%*%ginv(t(X)%*%solve(SB)%*%X)%*%t(X)%*%solve(SB)%*%ybar.B

    Q[j] <- ifelse(SIB>SI, 1, 0)
  }
  return(sum(Q)/length(Q))  ##p-value
}


#' multiple comparisons of the levels of AB interactions
#'
#'Using Parametric Bootstrap to simulate a distribution and find a p-value for the test
#'
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@usage Q.ABmc(L=5000, ns, means, s2, alpha=0.05, a, b, c)
#'@param alpha	significant level
#' @param ns	sample size for each group
#' @param means	sample mean for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return Q.crit:	The (1- alpha) percentile of the simulated distribution.
#' @return res.df: A dataframe containing the differences between each pair of factor levels, standard errors,
#' confidence interval for the differences, test statistic for each pair, p-value, and indicator of whether the
#' difference was statistically significant for each pair.
#' @return ybarij: estimated group mean for level i, j
#' @return var.YAB: estimated variance for each group mean
#' @return Q.test: largest test statistic from all pairs
#'
#'
#'@examples
#'
#'#'# We need elapsed time less than 5 seconds to publish
#'# the package with this example. Hence, there are # before
#'#a couple of PB algrithms. Remove # before you run the code.
#'
#'#Note that when running the example, the user should get similar p-values to the ones commented
#'# in the example, but they will not be identical.
#'attach(potato)
#'regime<-factor(regime)
#'variety<-factor(variety)
#'temp<-factor(temp)
#'
#'
#'#there are two levels for each factor, so a=b=c=2
#'library(Rmisc)
#'summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))
#'
#'
#'#need to extract group sizes (ns), group var's (s2), means (ybars) for function
#'pot.ns <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$N
#'pot.means <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$leak
#'pot.s2 <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$sd^2
#'#alg.ABC(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#'#0.1626, pvalue, so we do not reject, so we should drop the three way term.
#'
#'
#'
#'
#'#alg.BC(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#'#0.202, not significant, so the regime:temp interaction is not significant
#'#to check the other two-way interactions we need to reorder the data so that
#'#the 'BC' term is either regime:variety or temp:variety
#'
#'pot.ns.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp","regime","variety"))$N
#'pot.means.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp","regime","variety"))$leak
#'pot.s2.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp", "regime", "variety"))$sd^2
#'#alg.BC(ns=pot.ns.TRV, ybars=pot.means.TRV,s2=pot.s2.TRV, a=2, b=2, c=2, L=5000)
#'#p=0, reject H_0.  the regime:variety interaction is significant
#'
#'
#'pot.ns.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime", "temp","variety"))$N
#'pot.means.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime","temp","variety"))$leak
#'pot.s2.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime", "temp", "variety"))$sd^2
#'#alg.BC(ns=pot.ns.RTV, ybars=pot.means.RTV,s2=pot.s2.RTV, a=2, b=2, c=2, L=5000)
#'#p=0.0.3652, do not reject, the temp:variety interaction is not significant
#'
#'
#'##next we can see if we are able to drop the main effect 'temp',
#'##not involved with the regime:variety int.
#'##temp is factor 'A' in the TRV model above.
#'
#'##algorithm 4 tests factor C when only AB interaction is present.
#'## so we need the order that makes 'temp' factor C
#'#the way we originally ordered it, to test the ABC interaction
#'
#'#alg.C.AB(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#'#p-value is 0.002, so we cannot drop the temp term.  the final model is
#'
#'#y = variety + regime + temp + variety:regime
#'
#'
#'
#'Q.ABmc(ns=pot.ns, means=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#'
#'
#'#     95%
#'#2.832115
#'
#'#$res.df
#'#     groups differences std.errs      ci.lo      ci.hi test.stat      p sig
#'# 1 1 - 1 2   -1.803638 1.821063  -6.961098   3.353822 0.9904314 0.7676   -
#'# 1 1 - 2 1  -22.501936 2.895282 -30.701708 -14.302164 7.7719316 0.0000   *
#'# 1 1 - 2 2   -1.248118 1.615681  -5.823913   3.327677 0.7725027 0.8696   -
#'# 1 2 - 2 1  -20.698298 2.781589 -28.576079 -12.820517 7.4411765 0.0000   *
#'# 1 2 - 2 2    0.555520 1.401787  -3.414501   4.525541 0.3962943 0.9766   -
#'# 2 1 - 2 2   21.253818 2.651678  13.743962  28.763674 8.0152341 0.0000   *
#'
#'#$ybarij
#'#         [,1]     [,2]
#'#[1,]  4.91472 6.718358
#'#[2,] 27.41666 6.162838
#'
#'#$var.YAB
#'#        [,1]      [,2]
#'#[1,] 1.980845 1.3354254
#'#[2,] 6.401815 0.6295805
#'
#'#$Q.test
#'#[1] 8.015234
#'
#' @export
Q.ABmc <- function(L=5000, ns, means, s2, alpha=0.05, a, b, c){
  ##get the ns, means and s2 in an array so we can identify the indices
 # library(plyr)

  ns.ind <- arrange(expand.grid(A=1:a, B=1:b, C=1:c), A,B)
 # ns.ind <- arrange(expand.grid(A=1:a, B=1:b, C=1:c), A,B)

  n.grp <- array(0, c(a,b,c))  ##array does not fill entries in the desired order.
  s2.grp <- array(0, c(a,b,c))
  means.grp <- array(0, c(a,b,c))
  for(i in 1:a){
    for(j in 1:b){
      for(k in 1:c){
        n.grp[i,j,k] = ns[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
        s2.grp[i,j,k] = s2[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
        means.grp[i,j,k] = means[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
      }
    }
  }

  ##Calculate weights vk for actual test stat and the PB pivot variabl
  vk <- rep(0, c)
  for(k in 1:c){
    vk[k] <- sum(n.grp[,,k])
  }
  v.wt.k <- vk/sum(ns) ##the weights in order of the k index

  #calculate estimated means (using the weights) for each level of AB for the test statistic
  ybarij <- matrix(0, a, b)
  var.YAB <- matrix(0, a, b)

  for(i in 1:a){
    for(j in 1:b){
      ybarij[i,j] <- sum(v.wt.k*means.grp[i,j,])
      var.YAB[i,j] <- sum(v.wt.k^2 * s2.grp[i,j,]/n.grp[i,j,])
    }
  }

  ybarijVect <- as.vector(ybarij)
  var.YABvect <- as.vector(var.YAB)

  Qtest.mat <- matrix(0,a*b,a*b)
  #we just fill in upper triangular part
  for (r in 1: ((a*b) -1))
    for (s in (r+1):(a*b)){
      Qtest.mat[r,s]<- abs(ybarijVect[r] - ybarijVect[s])/sqrt(var.YABvect[r] + var.YABvect[s])
    }
  Q.test <- max(Qtest.mat)

  ##calculate the parts of the PB pivot variable
  Q <- rep(0, L)
  for(l in 1:L){  ##calculate the bootstrap means and sample variances
    y.B <- rep(0, length(means))
    s2.B <- rep(0, length(s2))

    for (j in 1:length(means)){
      y.B[j]<- rnorm(1, 0, sqrt(s2[j]/ns[j]))
      s2.B[j] <- rchisq(1, df=(ns[j]-1))*s2[j]/(ns[j]-1)
    }#end the j loop

    #put the bootstrap means and s2's in indexed arrays
    s2B.grp <- array(0, c(a,b,c))
    meansB.grp <- array(0, c(a,b,c))
    for(i in 1:a){
      for(j in 1:b){
        for(k in 1:c){
          s2B.grp[i,j,k] = s2.B[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
          meansB.grp[i,j,k] = y.B[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
          n.grp[i,j,k] = ns[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
        }
      }
    }

    #now Q will be the PB analogy of the Q.test above, use same weights
    yB.barij <- matrix(0, a, b)
    varB.YAB <- matrix(0, a, b)

    for(i in 1:a){
      for(j in 1:b){
        yB.barij[i,j] <- sum(v.wt.k*meansB.grp[i,j,])
        varB.YAB[i,j] <- sum(v.wt.k^2 * s2B.grp[i,j,]/n.grp[i,j,])
      }
    }

    yB.barijVect <- as.vector(yB.barij)
    varB.YABvect <- as.vector(varB.YAB)

    Qmat <- matrix(0,a*b,a*b)
    #we just fill in upper triangular part
    for (r in 1: ((a*b) -1))
      for (s in (r+1):(a*b)){
        Qmat[r,s]<- abs(yB.barijVect[r] - yB.barijVect[s])/sqrt(varB.YABvect[r] + varB.YABvect[s])
      }

    Q[l] <- max(Qmat)
  } #end l loop that has L reps

  Q.crit <-quantile(Q, 1-alpha)
  #  list(Q.crit = Q.crit, Q.test = Q.test)
  #modify function to return list of Q.crit, differences, SE and CI like the Tukey's test functions
  #ybarij is matrix of the est weighted means
  ybarij.v <- as.vector(t(ybarij)) #vector of weighted means in order y_11. , y_12. , ..., y_ab.
  varYAB.v <- as.vector(t(var.YAB)) #vector of variance est of each weighted mean in same order
  ab <- a*b
  diffs <- matrix(0, ab, ab)
  SEs <-  matrix(0, ab, ab)
  diff.inds <-  matrix(0, ab, ab)
  inds <- paste(expand.grid(1:a, 1:b)$Var2, expand.grid(1:a, 1:b)$Var1)
  for (k in 1:(ab-1)){
    for (i in 1:(ab-k)){
      diffs[i,i+k] <- ybarij.v[i] - ybarij.v[i+k]  #diff bw means (i, i+k)
      SEs[i, i+k] <- sqrt(varYAB.v[i] + varYAB.v[i+k])  #se corresponding to the diff's above
      diff.inds[i, i+k] <- paste(inds[i], "-", inds[i+k])
    } #end i loop
  } #end k loop
  diffs.v <- t(diffs)[lower.tri(t(diffs))]  #put them back in a list
  SEs.v <- t(SEs)[lower.tri(t(SEs))]
  diff.inds.v <- t(diff.inds)[lower.tri(t(diff.inds))]
  ci.lo <- diffs.v - SEs.v*Q.crit
  ci.hi <- diffs.v + SEs.v*Q.crit
  test.stat <- abs(diffs.v)/SEs.v
  signif <- ifelse(test.stat>Q.crit, "*", "-")
  pval<- 0
  for(i in 1:length(test.stat)){
    pval[i] <-	length(which(Q>test.stat[i]))/L
  }
  res.df <- data.frame(groups=diff.inds.v, differences=diffs.v, std.errs = SEs.v, ci.lo, ci.hi, test.stat=test.stat, p=pval, sig=signif)
  list(Q.crit=Q.crit, res.df=res.df, ybarij = ybarij, var.YAB=var.YAB, Q.test=Q.test)
}


#'for three-way ANOVA that has no significant interaction terms to test main effects
#'
#'Using Parametric Bootstrap to simulate a distribution and find a p-value for the test
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@usage alg.C(ns, ybars, s2, a, b, c, L)
#' @param ns	sample size for each group
#' @param ybars	sample mean for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return a simulated p-value for testing a main effect
#'
#'@examples
#'
#'#See Q.Amc
#'
#' @export
alg.C <- function(ns, ybars, s2, a, b, c, L){
  S <- diag(s2/ns) ##make S matrix
  ##make terms for X matrix
  J.abc <- rep(1, a*b*c)
  I.a <- diag(a)
  I.b <- diag(b)
  I.c <- diag(c)
  J.bc <- rep(1, b*c)
  J.a <- rep(1, a)
  J.b <- rep(1,b)
  J.c <- rep(1,c)
  I.ab <- diag(a*b)
  I.bc <- diag(b*c)
  X <- as.matrix(cbind( J.abc, kronecker(I.a, J.bc), kronecker(J.a, kronecker(I.b, J.c))))
  #test statistic

  SI <- t(ybars)%*%solve(S)%*%ybars -
    t(ybars)%*%solve(S)%*%X%*%ginv(t(X)%*%solve(S)%*%X)%*%t(X)%*%solve(S)%*%ybars
  ##Q, counts how many times test stat is less than PB pivot variable
  Q <- NULL
  for(j in 1:L) {
    ybar.B <- NULL
    S2B <-NULL
    for (i in 1:length(ybars)) {
      ybar.B[i] <- rnorm(1, mean=0, sd=sqrt(s2/ns)[i])  ##create bootstrap mean vector
      S2B[i] <- rchisq(1, df=(ns[i]-1)) * s2[i]/(ns[i]-1) ##create bootstrap variances vector
    }
    SB <- diag(S2B/ns)
    ##PB variable:
    SIB <-  t(ybar.B)%*%solve(SB)%*%ybar.B -
      t(ybar.B)%*%solve(SB)%*%X%*%ginv(t(X)%*%solve(SB)%*%X)%*%t(X)%*%solve(SB)%*%ybar.B
    Q[j] <- ifelse(SIB>SI, 1, 0)
  }
  return(sum(Q)/length(Q))  ##p-value
}

#'PB multiple comparisons of the levels of factor A (output like TukeyHSD)
#'
#'Using Parametric Bootstrap to simulate a distribution for the multiple comparisons and calculate a test stat
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@importFrom lmtest bptest
#'@usage Q.Amc(L=5000, ns, means, s2, alpha=0.05, a, b, c)
#'@param alpha	significant level
#' @param means	sample mean for each group
#' @param ns	sample size for each group
#' @param s2	sample variance for each group
#' @param a	Number of levels for factor A
#' @param b	Number of levels for factor B
#' @param c	Number of levels for factor C
#' @param L	Number of simulated values for the distribution
#' @return Q.crit: The (1-alpha) percentile of the simulated distribution.
#' @return Q.test: largest test statistic from all pairs
#' @return res.df: A dataframe containing the differences between each pair of factor levels, standard
#' errors, confidence interval for the differences, test statistic for each pair, p-value, and indicator
#' of whether the difference was statistically significant for each pair.
#'
#'@examples
#'# We need elapsed time less than 5 seconds to publish
#'# the package with this example. Hence, there are # before
#'#a couple of PB algrithms. Remove # before you run the code.
#'
#'#function to make everything but the response a factor
#'
#'make.factor <- function(dataset, fact.cols){
#'for( i in fact.cols){
#'  dataset[,i] <- factor(dataset[,i])
#'}
#'return(dataset)
#'}
#'
#'barley_ex <- make.factor(barleyh20, 1:5)
#'##this dataset has 4 factors, ignore year
#'
#'library(Rmisc)
#'
#'library(MASS)
#'
#'#library(lmtest)
#'
#'summarySE(barley_ex, "wt", c("genotype", "site", "time"))
#' #ignore year, note that the data are balanced
#'
#'summary(barley_ex$wt)
#'mod1 <- lm(wt~genotype*site*time, data=barley_ex)
#'anova(mod1)
#'plot(mod1$fit, mod1$resid)
#'qqnorm(mod1$resid)
#'shapiro.test(mod1$resid)
#'
#'boxcox(mod1, lambda=seq(-4, -2, by=0.1)) #lambda approx -3.5
#'mod2 <- lm(wt^(-3.5)~genotype*site*time, data=barley_ex)
#'plot(mod2$fit, mod2$resid) #worse?
#'
#'#go with untransformed data? drop 3way term
#'mod3 <- lm(wt~genotype + site + time + genotype:site + genotype:time + site:time, data=barley_ex)
#'anova(mod3) #site:time ns
#'
#'anova(lm(wt~genotype + site + time + genotype:site + genotype:time, data=barley_ex))
#'anova(lm(wt~genotype + site + time + genotype:site, data=barley_ex))
#'
#'anova(lm(wt~genotype + site + time, data=barley_ex))
#'anova(lm(wt~site + time, data=barley_ex))
#'anova(lm(wt~time, data=barley_ex))
#'TukeyHSD(aov(wt  ~ time, data=barley_ex))  #all sig except 35-30 and 20-15 (0.0569)
#'
#'
#'###use PB methods
#'summarySE(barley_ex, "wt", c("genotype", "site", "time"))
#'#note that the data are balanced
#'
#'#need to extract group sizes (ns), group var's (s2), means (ybars) for function
#'barley.ns <- summarySE(barley_ex, "wt", c("genotype", "site", "time"))$N
#'barley.means <- summarySE(barley_ex, "wt", c("genotype", "site", "time"))$wt
#'barley.s2 <- summarySE(barley_ex, "wt", c("genotype", "site", "time"))$sd^2
#'
#'#alg.ABC(ns=barley.ns, ybars=barley.means,s2=barley.s2, a=2, b=2, c=7, L=5000)
#'#p=0.9996, can drop 3way term
#'
#'#can we drop the site:time int term?
#'#alg.BC(ns=barley.ns, ybars=barley.means,s2=barley.s2, a=2, b=2, c=7, L=5000)
#'#p=0.9998, drop
#'
#'#reorder data to make the different two-way terms
#'barleyTSG.ns <- summarySE(barley_ex, "wt", c("time", "site", "genotype"))$N
#'barleyTSG.means <- summarySE(barley_ex, "wt", c("time", "site", "genotype"))$wt
#'barleyTSG.s2 <- summarySE(barley_ex, "wt", c("time", "site", "genotype"))$sd^2
#'
#'#alg.BC(ns=barleyTSG.ns, ybars=barleyTSG.means, s2=barleyTSG.s2, a=7, b=2, c=2, L=5000)
#'#p=0.9988, drop site:genotype
#'
#'#reorder to SGT, can we drop genotype:time?
#'barleySGT.ns <- summarySE(barley_ex, "wt", c("site", "genotype", "time"))$N
#'barleySGT.means <- summarySE(barley_ex, "wt", c("site", "genotype", "time"))$wt
#'barleySGT.s2 <- summarySE(barley_ex, "wt", c("site", "genotype", "time"))$sd^2
#'
#'#alg.BC(ns=barleySGT.ns, ybars=barleySGT.means,s2=barleySGT.s2, a=2, b=2, c=7, L=5000)
#'#p=0.9976, drop

#can we drop main effects?
#'#alg.C(ns=barley.ns, ybars=barley.means,s2=barley.s2, a=2, b=2, c=7, L=5000) #GST
#'#p=0, time has signif efffect, same conclusion as F-test
#'
#'#alg.C(ns=barleyTSG.ns, ybars=barleyTSG.means, s2=barleyTSG.s2, a=7, b=2, c=2, L=5000)
#'#p=0.9996 no signif effect of genotype
#'
#'
#'##site?
#'barleyGTS.ns <- summarySE(barley_ex, "wt", c("genotype", "time", "site"))$N
#'barleyGTS.means <- summarySE(barley_ex, "wt", c("genotype", "time", "site"))$wt
#'barleyGTS.s2 <- summarySE(barley_ex, "wt", c("genotype", "time", "site"))$sd^2
#'
#'#alg.C(ns=barleyGTS.ns, ybars=barleyGTS.means, s2=barleyGTS.s2, a=2, b=7, c=2, L=5000)
#'#p=0.9998, site is NS
#'
#'#multiple comparisons
#'#this function tests all pairwise comparisons of the levels of factor A,
#'# so we use the TSG order
#'Q.Amc(L=5000, ns=barleyTSG.ns, means=barleyTSG.means, s2=barleyTSG.s2,
#'alpha=0.05, a=7, b=2, c=2)
#'
#'#Demonstrate the method on unbalanced data by collapsing time into L, M, H
#'#barley_ex$time2 <- "M"
#'#barley_ex$time2 <- ifelse(as.numeric(barley_ex$time) <=2, "L", barley_ex$time2)
#'#barley_ex$time2 <- ifelse(as.numeric(barley_ex$time) >=6, "H", barley_ex$time2)
#'#barley_ex$time2 <- as.factor(barley_ex$time2)
#'
#'#still pretty balanced, separate the lowest level
#'#barley_ex$time2 <- ifelse(as.numeric(barley_ex$time) <2, "LL", barley_ex$time2)
#'#barley_ex$time2 <- as.factor(barley_ex$time2)
#'
#'#summarySE(barley_ex, "wt", c("genotype", "time2", "site"))
#'#two of the bigger groups N=12 have larger var
#'#anova(lm(wt ~genotype*time2*site, data=barley_ex))
#'#still looks like time2 the only sig factor
#'
#'#library(lmtest)
#'#bptest(lm(wt ~genotype*time2*site, data=barley_ex)) #violates
#'
#'#mod.un <- lm(wt ~genotype*time2*site, data=barley_ex)
#'#plot(mod.un$fit, mod.un$resid)
#'#qqnorm(mod.un$resid)
#'
#'#boxcox(lm(wt ~genotype*time2*site, data=barley_ex), lambda=seq(-6, -4, by=0.1))
#'#lambda = -4.5
#'
#'
#'#the above transformations didn't work so just try the untransformed data
#'#the three way interaction term was not significant
#'#anova(lm(wt ~genotype+ time2 + site + genotype:time2 + genotype:site + time2:site, data=barley_ex))
#'#anova(lm(wt ~genotype+ time2 + site + genotype:time2 + genotype:site, data=barley_ex))
#'#anova(lm(wt ~genotype+ time2 + site + genotype:time2, data=barley_ex))
#'#anova(lm(wt ~genotype+ time2 + site, data=barley_ex))
#'#anova(lm(wt ~genotype+ time2, data=barley_ex))
#'#anova(lm(wt ~time2, data=barley_ex))
#'
#'#TukeyHSD(aov(wt ~time2, data=barley_ex)) #all pairs significantly different
#'
#'
#'#PB methods
#'#barleyGST2.ns <- summarySE(barley_ex, "wt", c("genotype", "site", "time2"))$N
#'#barleyGST2.means <- summarySE(barley_ex, "wt", c("genotype", "site", "time2"))$wt
#'#barleyGST2.s2 <- summarySE(barley_ex, "wt", c("genotype", "site", "time2"))$sd^2
#'
#'
#'#alg.ABC(ns=barleyGST2.ns, ybars=barleyGST2.means,s2=barleyGST2.s2, a=2, b=2, c=4, L=5000)
#'#p=0.9734, can drop 3way term
#'
#'
#'#alg.BC(ns=barleyGST2.ns, ybars=barleyGST2.means,s2=barleyGST2.s2, a=2, b=2, c=4, L=5000)
#'#p=0.94, can drop site:time2
#'
#'#barleySGT2.ns <- summarySE(barley_ex, "wt", c("site","genotype",  "time2"))$N
#'#barleySGT2.means <- summarySE(barley_ex, "wt", c("site", "genotype", "time2"))$wt
#'#barleySGT2.s2 <- summarySE(barley_ex, "wt", c("site","genotype", "time2"))$sd^2
#'
#'#alg.BC(ns=barleySGT2.ns, ybars=barleySGT2.means,s2=barleySGT2.s2, a=2, b=2, c=4, L=5000)
#'#p=0.9952, can drop genotype:time2
#'
#'
#'#barleyTSG2.ns <- summarySE(barley_ex, "wt", c("time2", "site","genotype"))$N
#'#barleyTSG2.means <- summarySE(barley_ex, "wt", c("time2","site", "genotype"))$wt
#'#barleyTSG2.s2 <- summarySE(barley_ex, "wt", c("time2","site","genotype"))$sd^2
#'
#'
#'#alg.BC(ns=barleyTSG2.ns, ybars=barleyTSG2.means,s2=barleyTSG2.s2, a=4, b=2, c=2, L=5000)
#'#p=0.9556, can drop site:genotype
#'
#'#alg.C(ns=barleyGST2.ns, ybars=barleyGST2.means,s2=barleyGST2.s2, a=2, b=2, c=4, L=5000)
#'#p=0, time still has significant effect
#'
#'#alg.C(ns=barleyTSG2.ns, ybars=barleyTSG2.means,s2=barleyTSG2.s2, a=4, b=2, c=2, L=5000)
#'#p=0.9716, genotype is not significant
#'
#'
#'#barleyTGS2.ns <- summarySE(barley_ex, "wt", c("time2","genotype", "site"))$N
#'#barleyTGS2.means <- summarySE(barley_ex, "wt", c("time2","genotype", "site"))$wt
#'#barleyTGS2.s2 <- summarySE(barley_ex, "wt", c("time2","genotype", "site"))$sd^2
#'
#'
#'#alg.C(ns=barleyTGS2.ns, ybars=barleyTGS2.means,s2=barleyTGS2.s2, a=4, b=2, c=2, L=5000)
#'#p=0.9904, site is not significant
#'
#'
#'##mult comparisons of factor A so we put time first
#'#Q.Amc(L=5000, ns=barleyTSG2.ns, means=barleyTSG2.means, s2=barleyTSG2.s2,
#'# alpha=0.05, a=4, b=2, c=2)
#'#all sig, agrees with Tukey's test
#'
#'#TukeyHSD(aov(wt ~time2, data=barley_ex))
#'
#' @export
Q.Amc <- function(L=5000, ns, means, s2, alpha=0.05, a, b, c){
  ##Calculate weights for actual test stat and the PB pivot variable

  ns.ind <- arrange(expand.grid(A=1:a, B=1:b, C=1:c), A,B)
  n.grp <- array(0, c(a,b,c))  ##array does not fill entries in the desired order.
  for(i in 1:a){
    for(j in 1:b){
      for(k in 1:c)
        n.grp[i,j,k] = ns[which(ns.ind$A==i & ns.ind$B==j & ns.ind$C==k)]
    }
  }
  v.weight <- matrix(0, b, c)
  for(j in 1:b){
    for(k in 1:c){
      v.weight[j,k] <- sum(n.grp[,j,k])
    }
  }
  vjk <- as.vector(t(v.weight/sum(ns)))  ##the weights in order of the j,k index
  #calculate factor level estimated means (using the weights) for the test statistic
  ybari <- rep(0,a)
  ni <- rep(0,a)
  var.YA <- rep(0, a)
  ni[1] <- sum(ns[1:(b*c)])
  ybari[1] <- sum(vjk*means[1:(b*c)])
  var.YA[1] <- sum(vjk^2 * (s2/ns)[1:(b*c)])
  for(i in 2:a){
    ybari[i] <- sum(vjk*means[(b*c*(i-1)+1):(i*b*c)])
    ni[i] <- sum(ns[(b*c*(i-1)+1):(i*b*c)])
    var.YA[i] <- sum(vjk^2 * (s2/ns)[(b*c*(i-1)+1):(i*b*c)])
  }
  Qtest.mat <- matrix(0,a,a)
  diffs.mat <- matrix(0, a, a)
  diff.inds <- matrix(0, a, a)
  SEs.mat <- matrix(0, a, a)
  #we just fill in upper triangular part
  for (r in 1: (a -1))
    for (s in (r+1):(a)){
      Qtest.mat[r,s]<- abs(ybari[r] - ybari[s])/sqrt(var.YA[r] + var.YA[s])
      diffs.mat[r,s] <- ybari[s] - ybari[r]
      SEs.mat[r, s] <- sqrt(var.YA[r] + var.YA[s])
      diff.inds[r, s] <- paste(r, "-", s)
    }
  Q.test <- max(Qtest.mat)
  ##calculate the parts of the PB pivot variable
  Q <- rep(0, L)
  for(i in 1:L){  ##calculate the bootstrap means and sample variances
    y.B <- rep(0, length(means))
    s2.B <- rep(0, length(s2))

    for (j in 1:length(means)){
      y.B[j]<- rnorm(1, 0, sqrt(s2[j]/ns[j]))
      s2.B[j] <- rchisq(1, df=(ns[j]-1))*s2[j]/(ns[j]-1)
    }#end the j loop
    #now Q will be the PB analogy of the Q.test above. we use the same ni’s
    yB.bari <- rep(0,a)
    var.YBA <- rep(0, a)
    yB.bari[1] <- sum(vjk*y.B[1:(b*c)])
    var.YBA[1] <- sum(vjk^2 * (s2.B/ns)[1:(b*c)])
    for(m in 2:a){
      yB.bari[m] <- sum(vjk*y.B[(b*c*(m-1)+1):(m*b*c)])
      var.YBA[m] <- sum(vjk^2 * (s2.B/ns)[(b*c*(m-1)+1):(m*b*c)])
    } #end m loop
    Qmat <- matrix(0,a,a)
    #we just fill in upper triangular part
    for (r in 1: (a -1))
      for (s in (r+1):a){
        Qmat[r,s]<- abs(yB.bari[r] - yB.bari[s])/sqrt(var.YBA[r] + var.YBA[s])
      }
    Q[i] <-max(Qmat)
  } #end i loop that has L reps
  Q.crit <-quantile(Q, 1-alpha)
  #modify to return output like TukeyHSD function
  Qtest.vec <- as.vector(t(Qtest.mat)[lower.tri(t(Qtest.mat))])
  diffs.vec <- as.vector(t(diffs.mat)[lower.tri(t(diffs.mat))])
  SEs.vec <- as.vector(t(SEs.mat)[lower.tri(t(SEs.mat))])
  diff.inds.vec <- as.vector(t(diff.inds)[lower.tri(t(diff.inds))])
  ci.lo <- diffs.vec -SEs.vec*Q.crit
  ci.hi <- diffs.vec + SEs.vec*Q.crit
  signif <- ifelse(Qtest.vec>Q.crit, "*", "-")
  pval <- 0
  for(i in 1:length(Qtest.vec)){
    pval[i] <- length(which(Q>Qtest.vec[i]))/L
  }
  res.df <- data.frame(groups=diff.inds.vec, differences=diffs.vec, std.errs = SEs.vec, ci.lo, ci.hi, test.stat=Qtest.vec, p=pval, sig=signif)
  list(Q.crit = Q.crit, Q.test = Q.test, res.df=res.df)
}


#'PB multiple comparisons of the levels of factor A aginst the control group
#'
#'Using Parametric Bootstrap to simulate a distribution for the multiple comparisons of treatment
#'groups against a control
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@importFrom DescTools DunnettTest
#'@importFrom dplyr summarise
#'@usage dunnett.PB(L, ns, means, s2, alpha)
#'@param alpha	significant level
#' @param means	sample mean for each group
#' @param ns	sample size for each group
#' @param s2	sample variance for each group
#' @param L	Number of simulated values for the distribution
#' @return D.crit: The (1 - alpha) percentile of the simulated distribution
#' @return result: The differences, confidence intervals for the difference, and p-values for comparisons of each
#' factor level vs. the control.
#'
#'@examples
#'
#'#This one gets a different result between the PB method and the traditional Dunnett's test.
#'#Constant variance assumption appears violated on residual plots.
#'#The breusch pagan test shows close to violating (p=0.0596) while levene's test (using #median)
#'#does not show a violation.  Data is mildly unbalanced.
#'
#'#Traditional Dunnett's test says group 4-6 are different from "control", while the PB method
#'# only identifies group 5 and 6.
#'#Group 4 has a larger variance than the others. The pooled variance/MSE could be too small for
#'#this group and lead to arificially large test statistic.
#'#MSE is 0.0004274 and the sample variance of group 4 is 0.00169.
#'
#'#The authors of the paper do not claim significance, they report the means and state
#'#that there is a delineation between 30 and 40 feet.  This seems true when looking at the
#'#means, but the measurements at 40 feet do have a larger variance than those at other depths.
#'
#'#Practical interpretation:
#'#If your goal was to get the most iron rich water from as shallow depth as possible,
#'#knowing that the surface (control?) was not rich enough, and you decided to go 40 feet
#'# deep, you may still get water that didn't have enough iron content for your purpose.
#'#Suppose you wanted at least 0.1 content; based on the means you might use 40 feet, but
#'# from their data, the measurements were:
#'#0.098, 0.074, 0.154
#'#So 2/3 of these samples wouldn't contain enough iron for your purpose.
#'
#'
#'library(DescTools)  ##Dunnett's test
#'library(lmtest) #BP test for constant variance
#'library(dplyr) #data manipulation
#'library(MASS)
#'
#'fedata$depth <- factor(fedata$depth)
#'
#'femod <- lm(Y~depth, data=fedata)
#'plot(femod$fit, rstandard(femod),
#'     main="Fitted-Residual Plot, One-Way ANOVA Model", sub="Iron Data")
#'     #appears to violate equal variance assumption
#'
#'bptest(femod) #close to violation
#'
#'#what about normality?
#'qqnorm(rstandard(femod), main="Normal QQ-Plot, Standardized Residuals",
#'sub="One-Way ANOVA Model, Iron Data")
#'shapiro.test(rstandard(femod)) #does not violate
#'
#'fe.sums <-fedata %>% group_by(depth) %>% summarise(means=mean(Y),
#'vars = var(Y), sd=sd(Y), ns=length(depth))
#'fe.sums
#'#summarySE(fedata, "Y", "depth") from Rmisc pkg does same thing as fe.sums above
#'
#'pbd.fe <- dunnett.PB(L=5000, ns=fe.sums$ns, means=fe.sums$means,
#'                            s2=fe.sums$vars, alpha=0.05)
#'
#'pbd.fe$result
#'pbd.fe$D.crit
#'
#'##grp 5 and 6 sig diff from group 1
#'DunnettTest(Y~depth, data=fedata)
#'
#'##Dunnett's test also says group 4 is different.
#'##group 4 has a larger variance than the others.
#'##pooled variance/MSE could be too small for this group
#'##and lead to arificially large test statistic.
#'
#'anova(femod)
#'#MSE is 0.0004274
#'#sample variance of group 4 is 0.00169
#'
#'#Use traditional Dunnett's test on transformed data for comparisons
#'
#'#attempt log transformation
#'fedata$logY <-log(fedata$Y)
#'felogmod <- lm(logY~depth, data=fedata)
#'bptest(felogmod) #still violates
#'#LeveneTest(logY~depth, data=fedata)
#'plot(felogmod$fit, rstandard(felogmod))
#'shapiro.test(rstandard(felogmod)) #still for normality
#'
#'#square root transformation?
#'fedata$srY <- sqrt(fedata$Y)
#'fesrmod <- lm(srY~depth, data=fedata)
#'bptest(fesrmod) #still violates, worse
#'plot(fesrmod$fit, rstandard(fesrmod))
#'
#'boxcox(femod)
#'#lambda=-0.2
#'
#'fedata$bcY <- with(fedata, (Y^(-0.2) - 1)/-0.2)
#'febcmod <- lm(bcY~depth, data=fedata)
#'#gives same F-statistic as lm(Y^(-0.2) ~ depth, data=fedata),
#'
#'
#'bptest(febcmod)  #still violates
#'plot(febcmod$fit, rstandard(febcmod))
#'
#'##mg/L is a proportion, try arcsin
#'fedata$asY <- asin(sqrt(fedata$Y))
#'feasmod <-lm(asY~depth, data=fedata)
#'bptest(feasmod) #still close to violating
#'
#'plot(felogmod$fit, rstandard(felogmod), main="Log Transform.",
#'     xlab="Fitted Values", ylab="Standardized Residuals")
#'plot(febcmod$fit, rstandard(febcmod), main="Box-Cox Transform.",
#'     xlab="Fitted Values", ylab="Standardized Residuals")
#'plot(feasmod$fit, rstandard(feasmod), main="Arcsin(Sq. Root) Transform.",
#'     xlab="Fitted Values", ylab="Standardized Residuals")
#'
#'#P-values of BP test are similar for log and box-cox, plots look a little better
#'##log transform may be considered simpler, so try that
#'
#'anova(felogmod)
#'DunnettTest(logY~depth, data=fedata)  #still identifies 40 feet and above
#'
#'DunnettTest(bcY~depth, data=fedata) #still identifies 40 feet and above
#'shapiro.test(rstandard(febcmod)) # normality still ok
#'#W = 0.95535, p-value = 0.4556
#'
#'
#'
#' @export
dunnett.PB <- function(L, ns, means, s2, alpha){
  D <- rep(0, L)
  r <- length(ns) #number of groups
  pairs.data <- rep(0, r)
  diffs <- rep(0, r)
  SEs <- rep(0, r)
  for(j in 1:r){
    diffs[j] <- means[j]-means[1]
    SEs[j] <- sqrt( (s2[1]/ns[1]) + (s2[j]/ns[j]))
    pairs.data[j] <- abs(means[j]-means[1])/sqrt( (s2[1]/ns[1]) + (s2[j]/ns[j]))
    #the first 'pairs' will be 0
  }
  test.stat <- max(pairs.data)
  pairs <- rep(0, r)
  ##storage vector for the differences between group means for bootstrap data
  for(i in 1:L){
    y.B <- rep(0, r)
    s2.B <- rep(0, r)
    for (j in 1:r){
      y.B[j]<- rnorm(1)*sqrt(s2[j]/ns[j])
      s2.B[j] <- rchisq(1, df=(ns[j]-1))*s2[j]/(ns[j]-1)
      pairs[j] <- abs(y.B[1]-y.B[j])/sqrt( (s2.B[1]/ns[1]) + (s2.B[j]/ns[j]))
      #the first one will be 0
    }
    D[i]<- max(pairs)
  }
  pvals <- rep(0, r)
  for(j in 1:r){
    pvals[j] <- length(which(D>pairs.data[j]))/L
  }
  D.crit <- quantile(D, 1-alpha)
  list(result = data.frame(diffs=diffs, low.CI=diffs - D.crit*SEs,
                           upr.CI = diffs+D.crit*SEs, pvals=pvals),
       D.crit = D.crit)
}


#'PB test for main effect of one-way ANOVA
#'
#'Using Parametric Bootstrap to simulate a distribution for the main effect of one-way ANOVA
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@importFrom dplyr summarise
#'@usage alg.A1(ns, ybars, s2, a,L)
#'@param a	level of the factor
#' @param ybars	sample mean for each group
#' @param ns	sample size for each group
#' @param s2	sample variance for each group
#' @param L	Number of simulated values for the distribution
#' @return the simulated p-value
#'
#'@examples
#'
#'#see Q.Amc_oneway
#'
#'
#' @export
alg.A1 <- function (ns, ybars, s2, a, L) {
  S <- diag(s2/ns)
  X <- rep(1, a)
  SI <- t(ybars) %*% solve(S) %*% ybars - t(ybars) %*% solve(S) %*%
    X %*% ginv(t(X) %*% solve(S) %*% X) %*% t(X) %*% solve(S) %*%
    ybars
  Q <- NULL
  for (j in 1:L) {
    ybar.B <- NULL
    S2B <- NULL
    for (i in 1:length(ybars)) {
      ybar.B[i] <- rnorm(1, mean = 0, sd = sqrt(s2/ns)[i])  #PB means
      S2B[i] <- rchisq(1, df = (ns[i] - 1)) * s2[i]/(ns[i] - 1) #PB variance
    }
    SB <- diag(S2B/ns)
    SIB <- t(ybar.B) %*% solve(SB) %*% ybar.B - t(ybar.B) %*%
      solve(SB) %*% X %*% ginv(t(X) %*% solve(SB) %*% X) %*%
      t(X) %*% solve(SB) %*% ybar.B  #PB test stat, calculated L times
    Q[j] <- ifelse(SIB > SI, 1, 0)
  }
  return(sum(Q)/length(Q)) #return p-value
}

#'PB multiple comparisons of factor A in one-way ANOVA
#'
#'Using Parametric Bootstrap to simulate a distribution for multiple comparison in one-way ANOVA
#'
#'
#'@importFrom stats rnorm t.test glm pchisq quasipoisson uniroot quantile rchisq
#'@importFrom graphics abline lines par
#'@importFrom MASS ginv
#'@importFrom plyr arrange
#'@importFrom Rmisc summarySE
#'@importFrom dplyr summarise
#'@usage Q.Amc_oneway(L,ns, means, s2, alpha)
#'@param alpha	significant level
#' @param means	sample mean for each group
#' @param ns	sample size for each group
#' @param s2	sample variance for each group
#' @param L	Number of simulated values for the distribution
#' @return the simulated p-value
#' @return D.crit: The (1 - alpha) percentile of the simulated distribution
#' @return res.df: The differences, confidence intervals for the difference, and p-values for comparisons of each
#' two factor levels.
#'@examples
#'
#'library(pbANOVA)
#'
#'data(fedata)
#'
#'fedata$depth <- factor(fedata$depth)
#'
#'library(Rmisc)
#'summarySE(fedata, "Y", "depth")
#'
#'feNs <- summarySE(fedata, "Y", "depth")$N
#'feYs <- summarySE(fedata, "Y", "depth")$Y
#'fes2 <- (summarySE(fedata, "Y", "depth")$sd)^2
#'
#'anova(lm(Y~depth, data=fedata))  #F-test significant
#'#we saw in the dunnett's example that the equal variance assumption is violated
#'
#'library(MASS) #need MASS for ginv function for all the interaction and main effects algorithms
#'alg.A1(ns=feNs, ybars=feYs, s2=fes2, a=6, L=5000)
#'#p=0.0038
#'
#'#multiple comparisons
#'Q.Amc_oneway(L = 5000, ns=feNs, means=feYs, s2=fes2, alpha = 0.05)
#'
#'#compare to Tukey's test
#'TukeyHSD(aov(Y~depth, data=fedata))
#'
#'#results agree only for some levels.
#'
#'
#' @export
Q.Amc_oneway <- function (L = 5000, ns, means, s2, alpha = 0.05) {
  a <- length(means)
  Qtest.mat <- matrix(0, a, a)
  diffs.mat <- matrix(0, a, a)
  diff.inds <- matrix(0, a, a)
  SEs.mat <- matrix(0, a, a)
  for (r in 1:(a - 1)) for (s in (r + 1):(a)) {
    Qtest.mat[r, s] <- abs(means[r] - means[s])/sqrt((s2[r]/ns[r]) + (s2[s]/ns[s])) #test stat from real data
    diffs.mat[r, s] <- means[s] - means[r]     #entries for output with differences and CI's
    SEs.mat[r, s] <- sqrt((s2[r]/ns[r]) + (s2[s]/ns[s]))
    diff.inds[r, s] <- paste(r, "-", s)
  }
  Q.test <- max(Qtest.mat)  #dist of PB test stats
  Q <- rep(0, L)
  for (i in 1:L) {
    y.B <- rep(0, length(means))
    s2.B <- rep(0, length(s2))
    for (j in 1:length(means)) {
      y.B[j] <- rnorm(1, 0, sqrt(s2[j]/ns[j]))
      s2.B[j] <- rchisq(1, df = (ns[j] - 1)) * s2[j]/(ns[j] - 1)
    }
    Qmat <- matrix(0, a, a)
    for (r in 1:(a - 1)) for (s in (r + 1):a) {
      Qmat[r, s] <- abs(y.B[r] - y.B[s])/sqrt((s2.B[r]/ns[r]) + (s2.B[s]/ns[s]))
    }
    Q[i] <- max(Qmat)
  }
  Q.crit <- quantile(Q, 1 - alpha)
  Qtest.vec <- as.vector(t(Qtest.mat)[lower.tri(t(Qtest.mat))])
  diffs.vec <- as.vector(t(diffs.mat)[lower.tri(t(diffs.mat))])
  SEs.vec <- as.vector(t(SEs.mat)[lower.tri(t(SEs.mat))])
  diff.inds.vec <- as.vector(t(diff.inds)[lower.tri(t(diff.inds))])
  ci.lo <- diffs.vec - SEs.vec * Q.crit
  ci.hi <- diffs.vec + SEs.vec * Q.crit
  signif <- ifelse(Qtest.vec > Q.crit, "*", "-")
  pval <- 0
  for (i in 1:length(Qtest.vec)) {
    pval[i] <- length(which(Q > Qtest.vec[i]))/L
  }
  res.df <- data.frame(groups = diff.inds.vec, differences = diffs.vec,
                       std.errs = SEs.vec, ci.lo, ci.hi, test.stat = Qtest.vec,
                       p = pval, sig = signif)
  list(Q.crit = Q.crit, Q.test = Q.test, res.df = res.df)
}
