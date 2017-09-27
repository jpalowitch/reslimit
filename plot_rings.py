from graph_tool.all import *
import random
Gring = load_graph("Gring.graphml")
Gmesh = load_graph("Gmesh.graphml")
random.seed(12345)
pos = sfdp_layout(Gring)
pos2 = sfdp_layout(Gmesh)
for i in range(0, Gring.num_vertices()):
  pos2[i] = pos[i]
graph_draw(Gring, pos = pos, output = "Gring.png")
graph_draw(Gmesh, pos = pos2, output = "Gmesh.png")
