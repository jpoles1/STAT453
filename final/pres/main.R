library(BOIN)
library(ggplot2)
library(gridExtra)
library(plot3D)
f1 = function(i){(i["x"]+i["y"])}
f2 = function(i){(i["x"]+i["y"]**2)}
f3 = function(i){((i["x"])**2+(.2*i["y"]))}
f4 = function(i){(i["x"]+i["y"]**2+rnorm(1, mean = 0, sd = .2))}
f5 = function(i){((i["x"]-5)**3+(i["y"]-1)**3)}
plot_sim = function(fn, maxtox, name="fn", output=0){
  side = seq(1,4,.1)
  siderange = c(min(side), max(side))
  xy = as.vector(expand.grid(side, side))
  colnames(xy) = c("x", "y")
  z = apply(xy, 1, fn)
  #scale 0 to 1
  z = scale(z, center = min(z), scale = max(z)-min(z))*maxtox
  ds = cbind(xy, z)
  z = matrix(z, ncol=length(side), byrow=TRUE)
  #p = plot_ly(z=~z) %>% add_surface() %>% layout(title=paste('Plot of', name))
  #p = ggplot(data=ds, aes(x=x, y=y, z=z))+geom_contour()
  if(output){
    persp3D(side, side, z, xlab="Gemcitabine Lvl", ylab="MK-8776 Lvl", main = paste('Plot of', name), zlim=c(0, 1))
  }
  return(ds)
}
side = c(1:4)
scenario_sim = function(fn, maxtox, xside=side, yside=side, output=0){
  xside = scale(xside, center = min(xside), scale = max(xside)-min(xside))
  yside = scale(yside, center = min(yside), scale = max(yside)-min(yside))
  xy = as.vector(expand.grid(xside, yside))
  colnames(xy) = c("x", "y")
  z = apply(xy, 1, fn)
  #scale 0 to 1
  z = scale(z, center = min(z), scale = max(z)-min(z))*maxtox
  ds = cbind(xy, z)
  z = t(matrix(z, ncol=length(xside), byrow=TRUE))
  if(output){
    persp(xside, yside, z, xlab="Gemcitabine Lvl", ylab="MK-8776 Lvl", zlim=c(0, 1))
  }
  return(z)
}
sims = c(f1, f2, f3, f4, f5)
names = c("f1", "f2", "f3", "f4", "f5")
outcomes = c()
scenarios = c()
for(i in seq_along(sims)){
  sim = sims[[i]]
  sim_name = names[i]
  outcome = plot_sim(sim, 1, sim_name, output=0)
  scenario = scenario_sim(sim, 1)
  scenarios = append(scenarios, list(scenario))
  outcomes = rbind(outcomes, cbind(outcome, i))
}