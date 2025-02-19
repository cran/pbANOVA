% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{Q.ABmc}
\alias{Q.ABmc}
\title{multiple comparisons of the levels of AB interactions}
\usage{
Q.ABmc(L=5000, ns, means, s2, alpha=0.05, a, b, c)
}
\arguments{
\item{L}{Number of simulated values for the distribution}

\item{ns}{sample size for each group}

\item{means}{sample mean for each group}

\item{s2}{sample variance for each group}

\item{alpha}{significant level}

\item{a}{Number of levels for factor A}

\item{b}{Number of levels for factor B}

\item{c}{Number of levels for factor C}
}
\value{
Q.crit:	The (1- alpha) percentile of the simulated distribution.

res.df: A dataframe containing the differences between each pair of factor levels, standard errors,
confidence interval for the differences, test statistic for each pair, p-value, and indicator of whether the
difference was statistically significant for each pair.

ybarij: estimated group mean for level i, j

var.YAB: estimated variance for each group mean

Q.test: largest test statistic from all pairs
}
\description{
Using Parametric Bootstrap to simulate a distribution and find a p-value for the test
}
\examples{

#'# We need elapsed time less than 5 seconds to publish
# the package with this example. Hence, there are # before
#a couple of PB algrithms. Remove # before you run the code.

#Note that when running the example, the user should get similar p-values to the ones commented
# in the example, but they will not be identical.
attach(potato)
regime<-factor(regime)
variety<-factor(variety)
temp<-factor(temp)


#there are two levels for each factor, so a=b=c=2
library(Rmisc)
summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))


#need to extract group sizes (ns), group var's (s2), means (ybars) for function
pot.ns <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$N
pot.means <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$leak
pot.s2 <- summarySE(potato, measurevar="leak", groupvars=c("variety","regime","temp"))$sd^2
#alg.ABC(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#0.1626, pvalue, so we do not reject, so we should drop the three way term.




#alg.BC(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#0.202, not significant, so the regime:temp interaction is not significant
#to check the other two-way interactions we need to reorder the data so that
#the 'BC' term is either regime:variety or temp:variety

pot.ns.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp","regime","variety"))$N
pot.means.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp","regime","variety"))$leak
pot.s2.TRV <- summarySE(potato, measurevar="leak", groupvars=c("temp", "regime", "variety"))$sd^2
#alg.BC(ns=pot.ns.TRV, ybars=pot.means.TRV,s2=pot.s2.TRV, a=2, b=2, c=2, L=5000)
#p=0, reject H_0.  the regime:variety interaction is significant


pot.ns.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime", "temp","variety"))$N
pot.means.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime","temp","variety"))$leak
pot.s2.RTV <- summarySE(potato, measurevar="leak", groupvars=c("regime", "temp", "variety"))$sd^2
#alg.BC(ns=pot.ns.RTV, ybars=pot.means.RTV,s2=pot.s2.RTV, a=2, b=2, c=2, L=5000)
#p=0.0.3652, do not reject, the temp:variety interaction is not significant


##next we can see if we are able to drop the main effect 'temp',
##not involved with the regime:variety int.
##temp is factor 'A' in the TRV model above.

##algorithm 4 tests factor C when only AB interaction is present.
## so we need the order that makes 'temp' factor C
#the way we originally ordered it, to test the ABC interaction

#alg.C.AB(ns=pot.ns, ybars=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)
#p-value is 0.002, so we cannot drop the temp term.  the final model is

#y = variety + regime + temp + variety:regime



Q.ABmc(ns=pot.ns, means=pot.means,s2=pot.s2, a=2, b=2, c=2, L=5000)


#     95\%
#2.832115

#$res.df
#     groups differences std.errs      ci.lo      ci.hi test.stat      p sig
# 1 1 - 1 2   -1.803638 1.821063  -6.961098   3.353822 0.9904314 0.7676   -
# 1 1 - 2 1  -22.501936 2.895282 -30.701708 -14.302164 7.7719316 0.0000   *
# 1 1 - 2 2   -1.248118 1.615681  -5.823913   3.327677 0.7725027 0.8696   -
# 1 2 - 2 1  -20.698298 2.781589 -28.576079 -12.820517 7.4411765 0.0000   *
# 1 2 - 2 2    0.555520 1.401787  -3.414501   4.525541 0.3962943 0.9766   -
# 2 1 - 2 2   21.253818 2.651678  13.743962  28.763674 8.0152341 0.0000   *

#$ybarij
#         [,1]     [,2]
#[1,]  4.91472 6.718358
#[2,] 27.41666 6.162838

#$var.YAB
#        [,1]      [,2]
#[1,] 1.980845 1.3354254
#[2,] 6.401815 0.6295805

#$Q.test
#[1] 8.015234

}
