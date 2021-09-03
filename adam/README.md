# Subfolder for Adam

## Resources

Some resources I have found useful:

* https://www.andrewheiss.com/blog/2020/12/01/ipw-binary-continuous/
* https://livefreeordichotomize.com/2019/01/17/understanding-propensity-score-weighting/#how-do-we-incorporate-a-propensity-score-in-a-weight
* https://ehsanx.github.io/TMLEworkshop/

## Targetted MLE

From: https://ehsanx.github.io/TMLEworkshop/tmle.html#tmle-steps.
This looks a bit more involved

1) Transformation of continuous outcome variable
2) Predict from initial outcome modelling: G-computation
3) Predict from propensity score model
4) Estimate clever covariate H
5) Estimate fluctuation parameter
6) Update the initial outcome model prediction based on targeted adjustment of the initial predictions using the PS model
7) Find treatment effect estimate
8) Transform back the treatment effect estimate in the original outcome scale
9) Confidence interval estimation based on closed form formula
