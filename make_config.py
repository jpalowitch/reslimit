from graph_tool.all import *
import numpy
import sys
import os
b = numpy.concatenate([numpy.repeat(0, 25), numpy.repeat(1, 25)])
probs = numpy.array([[0.9, 0.1], [0.1, 0.9]])
probs = probs / sum(sum(probs)) * pow(150, 2)
out_degs = numpy.random.exponential(1, size = [1, 50])
g = graph_tool.generation.generate_sbm(b, probs, out_degs, out_degs, directed = False)
g_config = graph_tool.generation.random_rewire(g)

