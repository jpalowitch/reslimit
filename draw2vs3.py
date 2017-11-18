from graph_tool.all import *
import numpy
import sys
import os
g = load_graph("2vs3.gml")
pos = graph_tool.draw.arf_layout(g)
graph_tool.draw.graph_draw(g, pos=pos)
