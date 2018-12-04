# CandN

## The project contains files for Dubinkina, V., Fridman, Y., Pandey, P., & Maslov, S. (2018). Alternative stable states in a model of microbial community limited by multiple essential nutrients. arXiv:1810.04726. https://doi.org/10.1101/439547 ##

Contents of the folder:

- __Dynamic_Stability_L2/__ - contains code which performs dynamic perturbation analysis for 2Cx2Nx4S model

- __Marriage_Algorithm.ipynb__ - implementation of marriage algorithm for enumeration of allowed and uninvadable states for a given set of lambdas

- __MonteCarlo_Feasibility.ipynb__ - implementation of Monte-Carlo spanning of a flux space to get feasibility regions of states for a given set of lambdas and yields

- __Yield_Variation.ipynb__ - implementation of Monte-Carlo spanning of a flux space to get feasibility regions of states for a given set of lambdas and various yields drawn from a distribution

- __analyze_2Cx2N.R__ - code which performs analysis of Monte-Carlo spanning results for 2Cx2Nx4S model (for one set of yields and number of yields drawn from a distribution)

- __analyze_6Cx6N.R__ - code which performs analysis of Monte-Carlo spanning results for 6Cx6Nx36S model

- __functions.R__ - some customized functions for plotting

- __plots/__
...

Additional files with numerical simulation results used in the above scripts could be obtained from: [Named Link]LINK
