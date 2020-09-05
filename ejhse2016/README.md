# Dynamic hedging considering the degree of risk aversion

Code for paper in [HSE Economic Journal](https://ej.hse.ru/data/2016/03/25/1127994289/%D0%9B%D0%B0%D0%BA%D1%88%D0%B8%D0%BD%D0%B0.pdf).

This paper studies the problem of calculation the dynamic hedge ratio for the portfolio consisted of two assets. 
Commonly it’s solved assuming that the investor’s risk aversion is infinite. 
Then the optimal hedge coefficient is equal to ratio of covariance of the hedged and hedging assets to the variance of the latter. 
It’s natural to assume that the optimal hedge ratio also depends on the investor’s attitude to risk. 
In this paper this fact is implemented via maximization of the investor’s expected utility, which depends on the portfolio return and variance. 
Consequently if, for example, prices move upwards, the optimal hedge ratio is less than under the assumption of absolute risk aversion and vice versa. 
In the paper eight portfolios, consisted of Russian blue-chip stocks and futures, are estimated. 
Multivariate volatility models GO-GARCH and cop-GARCH are applied to estimate the conditional covariances and variances of hedged portfolio returns. 
There are additional parameters in the error term distribution, including skewness parameter, due to the existence of asymmetry effects in the financial assets returns’ distribution [Kroner, Ng, 1998].  
The hedge effectiveness is estimated on the out-of-sample period using the maximum attainable risk reduction, unconditional variance of hedged portfolio returns and financial result. 
It’s shown that in six cases cop-GARCH surpasses GO-GARCH in hedge. 
Including the degree of risk aversion in the investor’s utility function together with above-mentioned volatility models allows to obtain hedge effectiveness from 24% to 65%.
