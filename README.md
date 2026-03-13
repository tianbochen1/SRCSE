# SRCSE
Spatially regularized collective spectra estimation and clustering.

We propose spatially regularized collective spectral estimation (SRCSE) for jointly estimating spectral density functions (SDFs) from multiple series observed over irregular spatial locations, where each series is either a one-dimensional (1D) time series or a two-dimensional (2D) spatial field. SRCSE uses a common set of basis functions to capture similarities among the SDFs in a low-dimensional space, and estimates the basis coefficients by minimizing a Kullback--Leibler (KL) divergence criterion with two regularization terms. These regularization terms are introduced to enforce the smoothness of the estimated SDFs and to account for spatial dependence among the series. Simulations show that SRCSE achieves more accurate estimation of SDFs than competing methods and yields more spatially coherent clustering when the estimated SDFs are used for clustering. In particular, SRCSE produces clearer and smoother cluster boundaries with fewer isolated misclassified points while remaining robust under contamination. We further illustrate the practical utility of SRCSE through two real-data applications.

'fn.r': R code for functions.
'sim.r': reproduce the simulations.
'spanish.r': reproduce the results in Spanish weather data analysis.
