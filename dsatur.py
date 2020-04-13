from tbarray import TenseBinArray
from collections import namedtuple
from datetime import datetime
from copy import copy

	
__DEFAULT_BLOCKSIZE__ = 100


def __dsatur(
	edges, 
	blocksize1,
	blocksize2,
	mode,
	timeout,
	**kwargs   
):
	'''
	 DEGREE SATURATION ALGORITHM 
	
	 MODES:
	
	 mode = 0 -> run greedy, use rough assumptions for best_found and target
	 mode = 1 -> run branch-n-bound, use provided values for best_found and target
	
	
	 INPUT:
	  -  edges                := array of (int, int) containing graph's vertices
      -  blocksize1           := parameter for TenseBinArray containers for storing arrays of length
                                 related to vertices count
	  -  blocksize2           := parameter for TenseBinArray containers for storing arrays of length
                                 related to possible colors count
	  -  target               := if algorithm finds solution having this number of colors it stops
	  -  mode                 := 0 or 1, defines type of algorithm; 0 is for greedy, 1 is for branch-and-bound
	  -  timeout	          := int or None, limit for algorithm's work time in seconds; when exceeded best solution
	                             found so far will be returned or None
	  -  [kwargs]             := additional parameters required for branch-and-bound mode
	     -  best_found        := initial value for branch-and-bound mode which will be updated by algorithm
	     -  target            := target number of colors; if solution with this number of colors is found,
	                             algorithm stops; if target is set to 1 algorithm will continue working until
                                 optimal solution is found or time limit reached
								 
	
	 CONSTANTS:
	
	  -  vtx                  := set of graph's vertices
	  -  vtx_count            := number of graphs's vertices     
	  -  adj_vtx_matr         := matrix of size [vtx_count * vtx_count], each row is a TenseBinArray;
	  	                         if vertices X and Y are adjacent then adj_vtx_matr[X][Y] == 1 is True
	  -  target               := target value from input from branch-and-bound mode or 1 for greedy mode
	  -  timeout              := parameter from input
		
	
	 VARIABLES:
	
	  -  best_coloring        := None or array of int of length vtx_count; holds best solution found
	  -  best_found           := min number of colors in graph's coloring found so far 
	  -  vtx_to_cls           := array of int or None of lenght vtx_count; if color is assigned to vertex
	 		                     then vtx_to_cls[vertex] == color is True, otherwise vtx_to_cls[vertex] is None
	  -  cls_groups           := matrix of size [vtx_count * best_found]; each row is a TenseBinArray;
	 		                     if color X was assigned to vertex Y then cls_groups[X][Y] == 1 is True
	  -  adj_cls_matr         := matrix of size [best_found * vtx_count], each row is a TenseBinArray;
	 		                     if vertices X and Y are adjacent then adj_cls_matr[X][Y] == 1 is True 
	  -  adj_cls_count_matr   := matrix of size [best_found * vtx_count]; adj_cls_count_matr[X][Y] contains 
	 		                     number of neighbors of vertex X colored with color Y
	  -  adj_cls_count_all    := array of int of lenght vtx_count; for each vertex it holds number different 
	 		                     colors used to color its neighbors
	  -  stack                := stack used to run algorithm
	  -  stack_top            := index of top of stack
	
	
	 METHODS AND FUNCTIONS:
	
	  -  time_exceeded()      := if timeout parameter provided in input checks if out of time else; if parameter
	                             not provided always returns False
	  -  dstaur_score(V)      := calculates saturation degree and number of neighbors of vertex needed to choose
	                             next vertex to color, returns (int, int)
	  -  color_vertex(V,C)    := assigns color C to vertex V and updates all relevant data
	  -  uncolor_vertex(V,C)  := marks vertex V as uncolored and updates all relevant data
	'''
	
	#================= CONSTANTS AND VARIABLES ===============================
	
	vtx                = set(x for x, _ in edges) | set(y for _, y in edges)
	vtx_count          = len(vtx)
	target             = 1 if mode == 0 else kwargs["target"]
	adj_vtx_matr       = []
	
	best_coloring      = None
	best_found         = vtx_count if mode == 0 else kwargs["best_found"]
	vtx_to_cls         = [None] * vtx_count
	cls_groups         = [TenseBinArray([0] * vtx_count, blocksize1) for _ in range(best_found)]
	adj_cls_matr       = [TenseBinArray([0] * best_found, blocksize2) for _ in range(vtx_count)]
	adj_cls_count_matr = [[0] * best_found for _ in range(vtx_count)]
	adj_cls_count_all  = [0] * vtx_count
	adj_vtx_count      = []
	stack              = [None] * vtx_count
	stack_top          = 0
	
	# populate adj_vtx_matr and adj_vtx_count
	
	for vertex in vtx:
		
		array = [0] * vtx_count
		array[vertex] = 1
		adj_count = 0
		
		for e1, e2 in edges:
			if e1 == vertex:
				array[e2] = 1
				adj_count += 1
			elif e2 == vertex:
				array[e1] = 1
				adj_count += 1
		
		adj_vtx_matr.append(TenseBinArray(array, blocksize1))
		adj_vtx_count.append(adj_count)
		
		
	#=================== METHODS AND FUNCTIONS ===============================
	
	if timeout is None:
		time_exceeded = lambda : False
	else:
		start_time = datetime.now()
		time_exceeded = lambda : (datetime.now() - start_time).total_seconds() > timeout
	
	dsatur_score = lambda v: (adj_cls_count_all[v], adj_vtx_count[v])
	
	def adjacent_vertices(vertex):
		res = adjacent = filter(
			lambda idx: vertex != idx and adj_vtx_matr[vertex][idx] == 1,
			range(vtx_count)
		)
		return res
	
	def color_vertex(vertex, color):
	
		vtx_to_cls[vertex] = color
		cls_groups[color][vertex] = 1
	
		for idx in adjacent_vertices(vertex):	
			
			if adj_cls_matr[idx][color] == 0:
				adj_cls_matr[idx][color] = 1
				adj_cls_count_all[idx] += 1
			
			adj_cls_count_matr[idx][color] += 1
	
	def uncolor_vertex(vertex, color):
				
		vtx_to_cls[vertex] = None
		cls_groups[color][vertex] = 0
		
		for idx in adjacent_vertices(vertex):
			
			adj_cls_count_matr[idx][color] -= 1
			
			if adj_cls_count_matr[idx][color] == 0:
				adj_cls_matr[idx][color] = 0
				adj_cls_count_all[idx] -= 1
	
	
	#========================= ALGORITHM ===================================
		
	# find first vertex to process and put it on the top of the stack
	
	vertex = max(range(vtx_count), key=dsatur_score)
	stack[stack_top] = (vertex, None, 0)
	
		
	while stack_top != -1 and not time_exceeded():
	
		# take vertex from the top of the stack
		# if it is colored - undo relevant changes
	
		vertex, color, colors_used = stack[stack_top]
		
		if  color is not None:
			uncolor_vertex(vertex, color)
			start_color = color + 1
		else:
			start_color = 0
			
		# choose new color for this vertex
		# if number of colors used is less than best_found
		# then do necessary updates else go to next iteration
			
		possible_used_colors = filter(
			lambda idx: (adj_vtx_matr[vertex] & cls_groups[idx]).is_zero(),
			range(start_color, colors_used)
		)
		
		not_used_color = colors_used
		new_color = next(possible_used_colors, not_used_color)

		if new_color == not_used_color:
			colors_used += 1
			
		if colors_used >= best_found:
			stack_top -= 1
			continue
		else:
			stack[stack_top] = (vertex, new_color, colors_used)
			color_vertex(vertex, new_color)
		
		# choose next vertex to color and put it on the top of the stack
		# if no vertices left then new best solution is found
		# save it and go to next iteration if needed for b-n-b or exit for greedy
		
		non_colored_vertices = filter(
			lambda v: vtx_to_cls[v] is None,
			range(vtx_count)
		)
			
		next_vertex = max(non_colored_vertices, key=dsatur_score, default=None)
		
		if next_vertex is None:
			
			best_found = colors_used
			best_coloring = copy(vtx_to_cls)
			
			if mode == 0:
				break
			elif best_found <= target:
				break
			else:
				stack_top -= 1
			
		else:
			stack_top += 1
			stack[stack_top] = (next_vertex, None, colors_used)	
			
	return None if best_coloring is None else (best_found, best_coloring)
	
def color_graph(
	edges, 
	blocksize1 = __DEFAULT_BLOCKSIZE__,
	blocksize2 = __DEFAULT_BLOCKSIZE__,
	mode       = "greedy",
	improve    = None,
	timeout    = None
):

	if mode == "greedy":
	
		res = __dsatur(
			edges, 
			blocksize1, 
			blocksize2, 
			mode       = 0,
			timeout    = timeout
		)
		
	elif mode == "bnb":
		
		# get upper bound by running greedy algorithm
		greedy = __dsatur(
			edges, 
			blocksize1, 
			blocksize2, 
			mode       = 0,
			timeout    = timeout
		)
		
		# run branch-n-bound algorithm with constraint for upper bound
		if greedy is not None:
			
			bnb = __dsatur(
				edges, 
				blocksize1, 
				blocksize2, 
				mode       = 1,
				best_found = greedy[0],
				target     = 1 if improve is None else greedy[0] - improve,
				timeout    = timeout
			)
			
			res = greedy if bnb is None else bnb
			
		else:
			res = None
		
	else:
		raise ValueError("available modes are: greedy, bnb (branch-and-bound)")
		
	return res
	