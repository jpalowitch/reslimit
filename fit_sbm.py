import sys
import numpy
from graph_tool.all import *
command_args = sys.argv
fn = command_args[1]
g = load_graph(fn)
N = g.num_vertices()
state = minimize_blockmodel_dl(g)
blocks = state.get_blocks()
membership = [None] * N
for i in range(0, N):
  membership[i] = blocks[i]
numpy.savetxt(fname = command_args[2], X = membership)
