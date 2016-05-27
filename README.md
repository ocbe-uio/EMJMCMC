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

﻿**********************************Read Me******************************************

1. EMJMCMC 1.2.tar.gz (GNU zipped tar) file includes the most recent version of EMJMCMC package.

2. To install the package the user would need to run install.packages("https://github.com/aliaksah/EMJMCMC2016/files/270429/EMJMCMC_1.2.tar.gz", repos = NULL, type="source"). Notice that some dependencies may be required (see for details http://aliaksah.github.io/EMJMCMC2016/).

3. /R/the_mode_jumping_package2.r containes R OOP code for MJMCMC algorithm used in EMJMCMC package.

4. /paper/appedix.pdf contains proofs of the ergodicity of MJMCMC procedure and pseudo codes for MJMCMC and local combinatorial optimizers.

5. /examples/BAS/ archive contains the original BAS package that is addressed in the experiments along with EMJMCMC.

6. /examples/Simulated Data (Example 1)/ contains data (simcen-x.txt, simcen-y.txt) and code (mode_jumping_package_class_simulated_bas_data_1906.r, mode_jumping_package_class_simulated_bas_data_3211.r) for the first experiment. Notice that /examples/Simulated Data/BAS includes BAS based replications for the same experiment.

7. /examples/US Data/ contains U.S. Crime data (simcen-x1.txt, simcen-y1.txt) and code (mode_jumping_package_class_crime_bas_data_1909.r, mode_jumping_package_class_crime_bas_data_3237.r) for the second experiment. Notice that /examples/US Data/BAS includes BAS based replications for the same experiment.

8. /examples/Simulated Logistic Data With Multiple Modes (Example 3)/ contains simulated data (sim3-X.txt, sim3-Y.txt) and code (mode_jumping_package_class_example3_5000.r, mode_jumping_package_class_example3_10000.r) for the third experiment. Notice that /examples/Simulated Logistic Data With Multiple Modes (Example 3)/BAS includes BAS based replications for the same experiment.

9. /examples/Epigenetic Data/ contains Arabadopsis genetic and epigenetic data (epigen.txt) and code (epigenetic data poisson regression with a random effect 350.r, epigenetic data poisson regression with a random effect 500.r) for the fourth experiment. precalculated.csv contains the precalculated models for the Poisson regression with an AR(1) random effect for the addressed part of genome.

10. /examples/Protein Activity Data/ contains the protein activity data (proteincen.txt) and code (Protein activity data.r) for the fifth experiment. Notice that /examples/Protein Activity Data/BAS includes BAS based replications for the same experiment.

11. For additional details and updates see http://aliaksah.github.io/EMJMCMC2016/.

**********************************The End******************************************

