from utils import *
from polygon import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sa_vertices = [(4, 6), (6, 4), (5, 0), (3, 0), (2, 4)]
# obstacles = [[(3.2, 1), (3.2, 2), (4.2, 2), (4.2, 1)]]
plt.plot(list(map(lambda x: x[0], sa_vertices+sa_vertices[:1])), list(map(lambda x: x[1], sa_vertices+sa_vertices[:1])))
# fig, ax = plt.subplots()
# patches = []
# sa_polygon = Polygon(sa_vertices, True)
# patches.append(sa_polygon)
# ob_polygon = Polygon(obstacles[0], True)
# patches.append(ob_polygon)
# p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
# colors = 100*np.random.rand(len(patches))
# p.set_array(np.array(colors))

# ax.add_collection(p)
# sa_vertices = [(1, 0), (1, 4), (4, 4), (4, 0)]
sa = SurveyArea(sa_vertices, [])
# for e in sa.edges:
# 	print(e)
# print('eq: ', linearEqsToAMF(["y - 1.0x = 2.0"]))
sa.generateBCells()
plt.show()

# print(max(sa.vertices))

