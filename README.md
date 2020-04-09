# mr.ios

Insturments associated with other phenotypes than the exposure lead MR estimate being biased (horizontal pleiotropy). This package is designed to predict which instruments are likely to be problematic due to pleiotropy and downweight them to be used in MR. 

--------

## IOS
Suppose we estimate the effect of X on Y. X is instrumented by genetic variants, which are likely to be associated with other phenotypes. We can use external data to predict which instrument are the most pleiotropic to down weight them. 
 1. Make score for every SNP to see how likely it is to be pleiotropic.
 2. Perform MR accounting for pleiotropy using IOS (e.g. higher IOS = lower weight) 

* Two basic approaches to estimate IOS
	- ios1: A summary of the SNP r2 with the other traits (r2_gu)
	- ios2: A summary of the ratio of r2_gu / r2_gx, where r2_gx is the variance explained by the SNP on the exposure. Estimates the index of suspicion, whereupon SNPs which have a larger effect on a set of traits given their effect on the exposure are deemed more suspicious

---------

## Main question
 1. Does IOS correlate with Qj in model 1?
 2. Does Q reduce in model 2 compared to model 1?
 3. How different are casual estimates (bxy) in model 1 and model 2?
 4. How much does SE change (i.e. does power improve in model 2)?
