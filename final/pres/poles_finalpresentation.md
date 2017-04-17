Assessing the BOIN Algorithm using Bivariate Functions for Generative Toxicity Scenarios
========================================================
author: Jordan Poles
date: 4/18/17
autosize: false

Why Generate Toxicity Scenarios?
========================================================

<br>
### Herein, I use procedural generation to describe several different toxicity scenarios using simple bivariate equations. 

#### This allows greater creative freedom and repeatability.
<br>

```r
source("main.R")
example_bivariate = function(i){((i["x"]-13)**3+i["y"]-((i["x"]*i["y"])**3/i["y"]**2))}
```

========================================================
<center>
![plot of chunk unnamed-chunk-2](poles_finalpresentation-figure/unnamed-chunk-2-1.png)
</center>

What Tools Did I Use?
========================================================

- Created a set of functions in R which make it simple to compose a toxicity scenario using any bivariate equation (the results of which are automatically scaled from 0-1)
- Receive back both a plot of the toxicity scenario in 3D as well as a set of data points in both matrix and plot form.
- Levels for each of the two drugs in question may be set, such that functions can take an numeric array of ordered dosage levels (1-4) or quanitities (20mg, 40mg, 80mg); these values are input into the equation as fractions of the largest value.
- Experimenter can introduce stochasticity by incorporation of a random variable into the generative formula, as we show later.

Single MTD - Trial 1
========================================================


Single MTD - Trial 1
========================================================
![Single MTD - Trial 1](img/1-1.png)

Single MTD - Trial 1
========================================================
<br>
<p style="line-height: 150%">
For my first trial, I had envisioned a simple scenario for the toxicity matrix, but the software exhibited some strange behaviour. 
<br><br>
It only ended up exploring the dosage combinations for a single drug (Gemcitabine), without investigating any of the toxicity characteristics of the second. Thus, I do not believe it found a reasonable MTD.
</p>

Single MTD - Trial 2
========================================================
![Single MTD - Trial 2](img/1-2.png)


Single MTD - Trial 2
========================================================
<br>
<p style="line-height: 150%">
In my second trial, the search proceeded normally, and produced a more complete toxicity matrix as would be expected.
<br><br>
The algorithm was also able to find a reasonable choice for the MTD in this scenario. Thus it performed significantly stronger in this case as compared to the first trial.
</p>


Waterfall MTD - Trial 1
========================================================
![Waterfall MTD - Trial 1](img/2-1.png)

Waterfall MTD - Trial 1
========================================================
<br>
<p style="line-height: 150%">
In this trial, the algorithm performed a very thourough search of the toxicity matrix. 
<br><br>
It appears to have used the available subjects for each subtrial effectively. It seems that the chosen MTDs are very reasonable given the hypothetical trial scenario.
</p>

Waterfall MTD - Trial 2
========================================================
![Waterfall MTD - Trial 2](img/2-2.png)

Waterfall MTD - Trial 2
========================================================
<br>
<p style="line-height: 150%">
In this trial, the early-stopping mechanism of the design kicked-in. The remainder of the trials were not used, but the algorithm did not explore at least one notable dosage combination: (2,2). 
<br><br>
Nonetheless, the chosen MTD values are not unreasonable, though they may not permit maximum efficacy.
</p>
