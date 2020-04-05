# README

# TP 1

In the folder TP1, compile : `gcc tp1.c -O9 -o tp1`

## Exercise 5 - BFS

### Connected components

run "./tp1 graph.txt cc print" 

print = 0 : the program outputs the sizes of each connected components and the fraction of nodes in the largest connected component 	

print = 1 : the program outputs only the fraction of nodes in the largest connected component 

### Lower bound to the diameter

run "./tp1 graph.txt lower_bound"


## Exercise 6 - Triangles

run "./tp1 graph.txt triangle"

# TP 2 
## Exercise 1 
In this exercise we use the "Epinions social network opinion" data http://snap.stanford.edu/data/soc-Epinions1.html.

To generate the results run 



# TP 3

In the folder TP3, compile : `gcc tp3.c -O9 -o tp3`

## Exercise 1 - Simple benchmark  

run "./tp3 create_graphs" 

The files `graph_0.8_0.4.txt` (ratio = 2),  `graph_0.8_0.2.txt`(ratio = 4),  `graph_0.8_0.1.txt` (ratio = 8),  `graph_0.8_0.08.txt` (ratio = 10),  `graph_0.8_0.05.txt` (ratio = 16),  `graph_0.8_0.025.txt` (ratio = 32) will be generated.
## Exercise 2 - Label propagation

run `./tp3 graph.txt comm_graph.txt`

`comm_graph.txt` is the output file. 
On each line, there will be a node and its community

## Exercise 3 - Experimental evaluation
Requirements : numpy, matplotlib, networkx and scikit-learn

Run all the jupyter notebook cells of `tp3_visualization_and_evaluation.ipynb`
