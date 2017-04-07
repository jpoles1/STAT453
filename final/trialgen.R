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
  p = plot_ly(z=~z) %>% add_surface() %>% layout(title=paste('Plot of', name))
  #p = ggplot(data=ds, aes(x=x, y=y, z=z))+geom_contour()
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
  outcome = plot_sim(sim, sim_name, 1)
  scenario = scenario_sim(sim)
  scenarios = rbind(scenarios, cbind(scenario, i))
  outcomes = rbind(outcomes, cbind(outcome, i))
}
ggplot(data=outcomes, aes(x=x,y=y,z=z))+geom_contour()+facet_wrap(~i)
scenarios