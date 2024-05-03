# Raugraf

Pangenomic diversity metrics

# Local Complexity/Node Radius

```
raugraf local-complexity \
	--graph mygraph.gfa \        # GFA from PGGB (or at least passed thru smooothxg and odgi)
	--jumps 20 \                 # Traverse up to 20 nodes away from the focal node
	--node-bed out.bedGraph      # Bedgraph of local complexity per node for each path in the pangenome
```
