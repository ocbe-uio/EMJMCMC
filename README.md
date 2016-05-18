# EMJMCMC2016
A package with mode jumping MCMC for Bayesian variable selection and averaging within GLMM

Generalized linear mixed models (GLMM) are addressed for inference and
prediction in a wide range of different applications providing a powerful
scientific tool for the researchers and analysts coming from different
fields. In most of these fields more and more sources of data are becoming
available introducing a variety of hypothetical explanatory variables for
these models to be considered. Selection of an optimal combination of these
variables is thus becoming crucial in a Bayesian setting. The posterior
distribution of the models can be viewed as a relevant measure for the model
evidence, based on the observed data. The number of models to select from is
exponential in the number of candidate variables, moreover the search space
in this context is often extremely non-concave and has numerous local
extrema or statistically speaking modes. Hence efficient search algorithms
have to be adopted for evaluating the posterior distribution within a
reasonable amount of time. In this paper we introduce and implement
efficient mode jumping MCMC algorithms for calculating posterior
probabilities of the models for generalized linear models with a random
effect. Marginal likelihoods of models, given the specific choice of priors
and any choice of covariates, can be efficiently calculated using the
integrated nested Laplace approximations approach (INLA) for the class of
models addressed, however for some particular cases exact results are also
available. We further apply the suggested algorithm to some simulated data,
the famous U.S. crime data, gene expression data, and real epigenetic data
and compare its performance to some of the existing approaches like BAS, RS
or MC3.
