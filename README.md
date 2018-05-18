## bayesIndirect
The bayesIndirect R package examines the indirect effect of a SNP on the outcome through the mediator in a Bayesian framework with a spike and slab prior.

## Installation
```
install.packages("devtools") #The devtools package must be installed first
install.packages("coda") #The coda package must be installed first

devtools::install_github("SharonLutz/bayesIndirect")
```

## Example
For the given dataset, one can test if the SNP acts on forced explanatory volume (FEV) through the intermediate phenotype, cigarette smoking (smoke) after adjusting for age and gender. The code below runs this analysis.
```
library(bayesIndirect)
?bayesIndirect # For details on this function and how to choose input variables

data("dataB")
snp<-dataB[,"snp"]
smoke<-dataB[,"smoke"]
fev<-dataB[,"fev"]
age<-dataB[,"age"]
gender<-dataB[,"gender"]
bayesIndirect(x=snp,k=smoke,y=fev,z=cbind(age,gender),nCov=2)
```

## Output
For this analysis we have the following output and we can see that the snp acts on FEV through the intermediate phenotype cigarette smoking since 0 is contained in the 95 percent confidence bands for Beta. Please see the reference below for more details concerning this method.
```
1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:
            Mean       SD  Naive SE Time-series SE
gamma1  0.895653 0.305713 1.081e-03      1.209e-03
beta    0.071090 0.050192 1.775e-04      3.564e-04
alpha1 -0.001617 0.001238 4.378e-06      1.007e-05
alpha2  0.050219 0.065589 2.319e-04      3.989e-04
sig2    1.083435 0.048633 1.719e-04      1.710e-04

2. Quantiles for each variable:
            2.5%       25%       50%        75%     97.5%
gamma1  0.000000  1.000000  1.000000  1.0000000 1.0000000
beta   -0.009462  0.033432  0.072774  0.1066164 0.1678433
alpha1 -0.004052 -0.002464 -0.001612 -0.0007666 0.0007859
alpha2 -0.077065  0.005823  0.050013  0.0942904 0.1788785
sig2    0.992504  1.049987  1.081790  1.1154274 1.1828478
```

## Reference
**Lutz SM**, Hokanson JE, Sharma S, Weiss S, Raby B, Lange C. (2013) On the Integration of Expression Profiles in Genetic Association Studies: A Bayesian Approach to Determine the Path from Gene to Disease. *Open Journal of Genetics*. 3(3).
