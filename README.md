# local_boot

This repository present a new network bootstrap procedure, termed local bootstrap, to estimate standard errors of network statistics. 

Supporting functions related to this proposed methods are provided under the "functions" folders. The main function:

```
zhu_nb_boot <- function(A,quantile_n=0,B,returns = "boot",method = "own", distance = "zhu", 
kowning_u=NULL, induced_sampling=TRUE, weighted=FALSE, getT=NULL,...)
```

`A` is the observed adjacency matrix. `quantile_n` is the relative size of neighbor sets. Normally, you would only need to provide these two arugument to get a list of bootstrapping graphs. 

Additionally, you can provide your own distance function via `distance` and provide the true/estimated distance between vertices through `knowing_u`. If nothing provided, the algorithm will use a upper bound estiamtion as default distance. Weighted graphs are also supported, you just need to set `weighted` equal to `TRUE`.

Meanwhile, we also provide the R scripts to generate tables in related manuscript. 
