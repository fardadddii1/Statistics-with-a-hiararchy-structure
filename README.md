# Automatic statistics with a hiararchy structure
In this code, we perform an statistical analysis pipeline following a hierarchy of methods in which differnt tests are being applied. 
The code starts with Anderson darling (AD) Normality test. Based on the AD output, if data is normally distributedm it runs Barlett equivariance test, and if it is not normally distribute, it runs Levene's test. Then based on th results from normality test and equivariance, we apply either T-tes, Wilcoxon , or Welch's tests. 
