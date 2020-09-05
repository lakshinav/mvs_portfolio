# Hedging and Risk Aversion on Russian Stock Market: Strategies Based on MGARCH and MSV Models

Code for the paper in [Proceedings](http://ceur-ws.org/Vol-2018/#paper-10) of the The Second Workshop on Computer Modelling in Decision Making (CMDM), 2017.

The paper studies the problem of dynamic hedge ratio calculation for the portfolio consisted of two assets - stock and futures on that stock. Commonly it's solved assuming that the investor's risk aversion is infinite, but it's natural to assume that the optimal hedge ratio also depends on the investors attitude to risk. 
We apply investor's expected utility to implement risk aversion coefficient in the calculation of hedge ratio. In the paper seventeen portfolios, consisted of Russian blue-chip stocks and futures, are estimated. 
Multivariate volatility models GO-GARCH, copula-GARCH, asymmetric DCC and parsimonious stochastic volatility model are applied to estimate the conditional covariances of hedged portfolio returns. 
The hedge effectiveness is estimated on the out-of-sample period using the maximum attainable risk reduction, financial result and investor's utility.
It's shown that for 60% portfolios ADCC surpasses the other models in hedging. 
Including the degree of risk aversion in the investor's utility function together with above-mentioned volatility models allows to reach hedge effectiveness of 88%.
