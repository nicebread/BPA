## ----setup, include=FALSE-----------------------------------------------------
	library(tint)
	knitr::opts_chunk$set(cache=TRUE, warnings = FALSE, messages=FALSE)
	# load all functions
	devtools::load_all()

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("nicebread/BFDA", subdir="package")

## ----eval=FALSE, fig.width=7, fig.height=4------------------------------------
#  sim.H1 <- BFDA.sim(expected.ES = 0.5, ...)
#  sim.H0 <- BFDA.sim(expected.ES = 0, ...)
#  
#  BFDA.analyze(sim.H1)
#  BFDA.analyze(sim.H0)
#  
#  plot(sim.H1)
#  plot(sim.H0)
#  
#  SSD(sim.H1)
#  SSD(sim.H0)

## ----sim, cache=TRUE, results='hide', message=FALSE---------------------------
sim.H1 <- BFDA.sim(expected.ES=0.5, type="t.between",
                   prior=list("Cauchy",list(prior.location=0, prior.scale=sqrt(2)/2)),
                   n.min=20, n.max=300, alternative="greater", boundary=Inf, B=1000,
                   verbose=TRUE, cores=1, stepsize = 10)

sim.H0 <- BFDA.sim(expected.ES=0, type="t.between",
                   prior=list("Cauchy", list(prior.location=0, prior.scale=sqrt(2)/2)),
                   n.min=20, n.max=300, alternative="greater", boundary=Inf, B=1000,
                   verbose=TRUE, cores=1, stepsize = 10)

## ----analyze------------------------------------------------------------------
BFDA.analyze(sim.H1, design="fixed", n=50, boundary=6)
BFDA.analyze(sim.H0, design="fixed", n=50, boundary=6)

## ----analyze2-----------------------------------------------------------------
BFDA.analyze(sim.H1, design="sequential", n.min=20, n.max=300, boundary=10)

## ----analyze3-----------------------------------------------------------------
BFDA.analyze(sim.H1, design="sequential", n.min=20, n.max=100, boundary=10)

## ----evDens, warning=FALSE, fig.width=7, fig.height=4-------------------------
evDens(BFDA.H1=sim.H1, BFDA.H0=sim.H0, n=20, boundary=c(1/6, 6), xlim=c(1/11, 31))

## ----SBF1, warning=FALSE, fig.width=7, fig.height=4---------------------------
plot(sim.H1, n.min=20, boundary=c(1/6, 6), n.trajectories = 60)

## ----SBF0, warning=FALSE, fig.width=7, fig.height=4---------------------------
plot(sim.H0, n.min=20, boundary=c(1/3, 3))

## ----SBF_nmax1, warning=FALSE, fig.width=7, fig.height=4----------------------
plot(sim.H1, n.min=20, n.max=80, boundary=c(1/5, 10))

## ----SBF_nmax0, warning=FALSE, fig.width=7, fig.height=4----------------------
plot(sim.H0, n.min=20, n.max=80, boundary=c(1/5, 10), forH1 = FALSE)

## ----SSD1, warning=FALSE------------------------------------------------------
SSD(sim.H1, power=.80, boundary=c(1/10, 10))

## ----SSD2, warning=FALSE------------------------------------------------------
SSD(sim.H0, alpha=.02, boundary=c(1/6, 6))

## ----paired.t, warning=FALSE, fig.width=7, fig.height=4-----------------------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
s1 <- BFDA.sim(expected.ES=0.4,
               prior=list("t", list(prior.location=0, prior.scale=sqrt(2)/2, prior.df=1)),
               n.min=50, stepsize=5, n.max=300, type="t.paired", design="sequential",
               alternative="greater", B=1000, cores=1, verbose=FALSE)
s0 <- BFDA.sim(expected.ES=0,
               prior=list("t", list(prior.location=0, prior.scale=sqrt(2)/2, prior.df=1)),
               n.min=50, stepsize=5, n.max=300, type="t.paired", design="sequential",
               alternative="greater", B=1000, cores=1, verbose=FALSE)

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(s1, design="sequential", boundary=10)
BFDA.analyze(s0, design="sequential", boundary=10)

BFDA.analyze(s1, design="sequential", boundary=6)
BFDA.analyze(s0, design="sequential", boundary=6)

plot(s1)

## ----corr_kappa_plot1, cache=TRUE, results='hide', fig.width=7, fig.height=7, echo=FALSE----
rho <- seq(-1, 1, by=.01)

par(mfrow=c(2, 3))
plot(rho, BFDA:::.priorRho(rho, kappa=1), main=bquote(kappa*" = 1"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRho(rho, kappa=1.5), main=bquote(kappa*" = 1.5"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRho(rho, kappa=2), main=bquote(kappa*" = 2"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRho(rho, kappa=3), main=bquote(kappa*" = 3"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRho(rho, kappa=0.75), main=bquote(kappa*" = 0.75"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRho(rho, kappa=0.001), main=bquote(kappa*" = 0.001"),
     ylim=c(0, 2.4), type="l", xlab="rho", ylab="Plausibility")

## ----corr_kappa_plot2, cache=TRUE, results='hide', fig.width=7, fig.height=7, echo=FALSE----
rho <- seq(0, 1, by=.01)

par(mfrow=c(2, 3))
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=1), main=bquote(kappa*" = 1"), 
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=1.5), main=bquote(kappa*" = 1.5"), 
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=2), main=bquote(kappa*" = 2"), 
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=3), main=bquote(kappa*" = 3"), 
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=0.75), main=bquote(kappa*" = 0.75"), 
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")
plot(rho, BFDA:::.priorRhoPlus(rho, kappa=0.001), main=bquote(kappa*" = 0.001"),
     ylim=c(0, 4.2), type="l", xlab="rho", ylab="Plausibility")

## ----corr_walkthrough, warning=FALSE, fig.width=7, fig.height=4---------------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
c1 <- BFDA.sim(expected.ES=0.21, prior=list("stretchedbeta",list(prior.kappa=2)),
               n.min=50, stepsize=10, n.max=300, B=1000, type="correlation",
               design="sequential", alternative="greater", cores=1, verbose=FALSE)
c0 <- BFDA.sim(expected.ES=0, prior=list("stretchedbeta", list(prior.kappa=2)),
               n.min=50, stepsize=10, n.max=300, B=1000, type="correlation",
               design="sequential", alternative="greater", cores=1, verbose=FALSE)

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(c1, design="sequential", boundary=10)
BFDA.analyze(c0, design="sequential", boundary=10)

BFDA.analyze(c1, design="sequential", boundary=6)
BFDA.analyze(c0, design="sequential", boundary=6)

plot(c1, boundary=c(1/10, 20), n.max=150)

## ----abtest_walkthrough, warning=FALSE, fig.width=7, fig.height=4-------------
#devtools::install_github("nicebread/BFDA", subdir="package")
#library(BFDA)

# do a sequential design analysis
ab1 <- BFDA.sim(expected.ES=2, prior=list("normal", list(prior.mean=0.5, prior.variance = 1)),
                n.min=50, stepsize=10, n.max=300, B=1000, type="abtest",
                design="sequential", alternative="two.sided", cores=1, verbose=FALSE,
                options.sample = list(effecttype = "OR"))
ab0 <- BFDA.sim(expected.ES=0, prior=list("normal", list(prior.mean=0.5, prior.variance = 1)),
                n.min=50, stepsize=10, n.max=300, B=1000, type="abtest",
                design="sequential", alternative="two.sided", cores=1, verbose=FALSE,
                options.sample = list(effecttype = "OR"))

# if no n.min and n.max is provided in the `BFDA.analyze` function,
# the values from the simulation are taken
BFDA.analyze(ab1, design="sequential", boundary=10)
BFDA.analyze(ab0, design="sequential", boundary=10)

BFDA.analyze(ab1, design="sequential", boundary=6)
BFDA.analyze(ab0, design="sequential", boundary=6)

plot(ab1, boundary=c(1/10, 20), n.max=150)

## ----grant--------------------------------------------------------------------
# We use the simulation from above. 
# Check the expected sample sizes for an evidential boundary of 10
a1 <- BFDA.analyze(sim.H1, design="sequential", n.min=20, boundary=10)

# --> see 80% quantile in output
a1

# Alternative approach: access stopping-ns directly
n_q80 <- ceiling(quantile(a1$endpoint.n, prob=.80))
n_q80

## ----grant2-------------------------------------------------------------------
a2.H1 <- BFDA.analyze(sim.H1, design="sequential", n.min=20, n.max=n_q80, boundary=10)
a2.H0 <- BFDA.analyze(sim.H0, design="sequential", n.min=20, n.max=n_q80, boundary=10)
a2.H1
a2.H0

