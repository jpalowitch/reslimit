from graph_tool.all import *
import numpy
import sys
import os
tmp_dir = "/tmp/RtmpBNNCFx" # sys.argv[1]
b = numpy.loadtxt(os.path.join(tmp_dir, "b.dat"))
probs = numpy.loadtxt(os.path.join(tmp_dir, "probs.dat"))
out_degs = numpy.loadtxt(os.path.join(tmp_dir, "out_degs.dat"))
g = graph_tool.generation.generate_sbm(b, probs, out_degs, directed = False)
