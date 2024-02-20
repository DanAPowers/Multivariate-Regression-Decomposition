### Multivariate-Regression-Decomposition
This module decomposes an observed group difference in means/proportions/counts/rates estimated from multivariate regression models into components reflecting group differences in characteristics and group differences in the effects of those characteristics. The methodology is described in D.Powers, H.Yoshioka and M-S. Yun -- mvdcmp: Multivariate Decomposition for Nonlinear Response Models, Stata Journal (2011). 11, Number 4, pp. 556–576.

Multivariate decomposition is widely used in research to quantify the contributions to group differences in average predictions from multivariate models. The technique utilizes the output from regression models to parcel out components of a group difference in a statistic, such as a mean or proportion, which can be attributed to compositional differences between groups (i.e., differences in characteristics or endowments) and to differences in the effects of characteristics (i.e., differences in the returns, coefficients, or behavioral responses). These techniques are equally applicable for partitioning change over time into components attributable to changing effects and changing composition.

This module implements regression decomposition for a variety of models. Data must be setup as a complex survey design (although this also accommodates simple designs). The svyglm procedure in T. Lumley's survey package is the workhorse of the procedure. The svyglm.nb procedure in D. Ludecke's sjstats package is used for negative binomial models. A robust estimator of the standard errors is the default for all models. A typical call is:

    call: decomp.model(formula, Asub, Bsub, scale=1, printit=FALSE, reverse=FALSE)
    
where Asub, and Bsub are subsets representing the two groups/time periods whose difference in 1st moments is to be decomposed. Options allows for scaling of the differences, outputing results to the terminal, and performing a reverse decomposition that swaps groups.

The following models are supported:

    decomp.linear      (linear regression decomposition)
    decomp.logit       (logit regression decomposition)
    decomp.probit      (probit regression decomposition)
    decomp.cloglog     (complementary log-log regression decomposition)
    decomp.poisson     (poisson regression decomposition)
    decomp.negbin      (negative binomial regression decomposition)

The example file mvdecompExample.R provides more detail using data from the NHANES.
