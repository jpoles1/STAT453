require("BOIN")
require("ggplot2")
n_sim = 5
shrink_matrix = function(x, min=0, max=1){
  
}
f1 = function(x){(x[1]+x[2])/8}
xy = as.vector(expand.grid(c(1:4), c(1:4)))
colnames(xy) = NULL
z = apply(xy, 1, f1)
ggplot(aes(x=c(1:4), y=c(1:4), z=z))
for(i in n_sim){
  p.true = matrix(c(0.02, 0.04, 0.08, 0.14, 0.08, 0.25, 0.42, 0.48, 0.25, 0.45,
0.50, 0.60), ncol=4, byrow=TRUE)
}
