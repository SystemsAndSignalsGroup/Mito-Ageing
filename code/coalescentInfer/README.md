# Inference Steps

The scripts in this folder will precalculate a Monte Carlo integration required to grid search over the likelihood given by the heirarchical model, then use this precalculated matrix to infer the posterior given the data provided. (Warning precalculating takes a long time, sample size can be reduced sacrificing the precision of the integration for speed).

1) Run the monteIntegralPrecalcTimeSplit.sh to generate folders of precalculated integrations.
2) Run concatIntegration.py to create a single file (logLikelySingleCell90.pkl) with all precalculated integrations in.
3) Run finalHeirarchy.sh to calculate the required marginal posteriors for plotting.
