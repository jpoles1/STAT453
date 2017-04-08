#STAT 453 - Exam Notes
###By: Jordan Poles
http://odin.mdacc.tmc.edu/~yyuan/index_code.html
##Question 1
_Use a stoping size equal to cohort size!_

_Create your own  toxicity scenarios._
####Part 1

Here, I simulate the operating characteristics of a BOIN trial design given 5 different meaningful toxicity scenarios of my own design. I decided to take a slightly novel approach, and use procedural generation to describe several different toxicity scenarios using simple multivariate equations.

I've created a set of functions in R which make it simple to functionally compose a toxicity scenario using any bivariate equation (the results of which are automatically scaled from 0-.8; assuming we have chosen a maximum dose that should not produce 100% toxicity), and receive back both a plot of the toxicity scenario in 3D as well as a set of data points (a toxicity level for each level of treatment) which may be input directly into the BOIN R/GUI software in order to discretely simulate the scenario.   

```{r, echo=FALSE}
require(BOIN)
require(ggplot2)
require(gridExtra)
require(plotly)
f1 = function(i){(i["x"]+i["y"])}
f2 = function(i){(i["x"]+i["y"]**2)}
f3 = function(i){((i["x"]-3)**2-i["y"])}
#f4 = function(i){((i["x"])**2+.1*i["y"])}
f4 = function(i){((i["x"]+3)**4*(i["y"]+3)**4)}
f5 = function(i){((i["x"]-3)**2+(i["y"]-3)**2)}
plot_sim = function(fn, name="fn", output=0){
  side = seq(1,4,.1)
  #side = c(1:4)
  siderange = c(min(side), max(side))
  xy = as.vector(expand.grid(side, side))
  colnames(xy) = c("x", "y")
  z = apply(xy, 1, fn)
  #scale 0 to 1
  z = scale(z, center = min(z), scale = max(z)-min(z))
  ds = cbind(xy, z)
  z = matrix(z, ncol=length(side), byrow=TRUE)
  #p = plot_ly(z=~z, type="surface") %>% layout(title=paste('Plot of', name))
  p = ggplot(data=ds, aes(x=x, y=y, z=z))+geom_contour()
  if(output){print(p)}
  return(ds)
}
scenario_sim = function(fn){
  side = c(1:4)
  siderange = c(min(side), max(side))
  xy = as.vector(expand.grid(side, side))
  colnames(xy) = c("x", "y")
  z = apply(xy, 1, fn)
  #scale 0 to 1
  z = scale(z, center = min(z), scale = max(z)-min(z))
  ds = cbind(xy, z)
  z = matrix(z, ncol=length(side), byrow=TRUE)
  return(z)
}
sims = c(f1, f2, f3, f4, f5)
names = c("f1", "f2", "f3", "f4", "f5")
outcomes = c()
scenarios = c()
for(i in seq_along(sims)){
  sim = sims[[i]]
  sim_name = names[i]
  outcome = plot_sim(sim, sim_name, 0)
  scenario = scenario_sim(sim)
  scenarios = rbind(scenarios, cbind(scenario, i))
  outcomes = rbind(outcomes, cbind(outcome, i))
}
#s = plot_sim(f1, output=1)
ggplot(data=outcomes, aes(x=x,y=y,z=z, colour = ..level..))+geom_contour()+facet_wrap(~i)
#scenarios
```

A summary of our different toxicity scenarios:
1) A linear relationship between both drugs. When either drug is at maximum dosage, the toxicity level is .5,
2)
3)
4)
5)

####Part 2

Here we create a realistic phase 1 clinical trial scenario and

_We want to know if there are any dose transitions recommended by the software which look a bit fishy?_
##Question 2
_Set target toxicity probability (under target probability) to .25_

_Use the sample size stuff that's in there by default_

##Question 3
- You should be able to assign a drug quanity for each dosage level.
- You should be able to import/export curves as a CSV list.
