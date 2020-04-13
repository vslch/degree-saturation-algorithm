# degree-saturation-algorithm
Implementation of degree saturation algorithm for graph coloring problem. Two modes are available - greedy and branch-and-bound. Algorithm involves sets intersection operations done via bitwise-logic with custom class `TenseBinArray`, see **tbarray.py**. 

# Usage

`color_graph` function in **dsatur.py**

## Parameters

- `edges` - list of tuples (int, int); represent graph's edges; 

- `blocksize1` - (optional) int, default is 100; parameter for storing binary data associated with vertices, for more  details see implementation;

- `blocksize2` - (optional) int, default is 100; parameter for storing binary data associated with color, for more details see implementations;

- `mode` - (optional) string, available modes are `"greedy"` and `"bnb"`, default is `"greedy"`; defines version of algorithm to run

- `timeout` - (optional) int or None, default is None; defines limit for algorithm's work time in seconds, when exceeded best solution found so far will be returned or None; 

- `improve` - (optional) int or None, default is None; optional parameter for `bnb` mode, if passed to the function then algorithm will stop working when solution with less or equal colors then `G - improve` is found where `G` stands for solution found by greedy approach; if solution with this property is not found, best found solution will be returned.  

## Returns

- tuple (num, color) where `num` is number of colors used and `color` is an array of ints so that `color[i]` is the color assigned to i'th vertex

## Example
```console
>>> from dsatur import color_graph
>>> color_graph([(4,3),(0,1),(1,2),(1,3)])
(2, [1, 0, 1, 1, 0])
```
