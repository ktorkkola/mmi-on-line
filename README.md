mmi-on-line
===========

Dimension reduction by mutual information maximization; the simple on-line learning version.

See appendix A of 
[Feature Extraction by Non-Parametric Mutual Information Maximization] (http://www.jmlr.org/papers/v3/torkkola03a.html)
Kari Torkkola; Journal of Machine Learning Research 3(Mar):1415-1438, 2003.


File pfMex.* (and mex-compiling it to your platform) is unnecessary unless you want to plot the information forces. 
Start by running script demo.m.
There are three important parameters:
* kernel width par.sigma (also try changing from Gaussian to a more heavy-tailed kernel for fun)
* learning rate decays linearly from par.initialEta to par.Finaleta
* par.batchSize how many random pairs of points sampled for each step


