from utils import *
from polygon import *

# sa_vertices = [(3, 0), (0, 3), (4, 6), (6, 4), (5, 0)]
sa_vertices = [(1, 0), (1, 4), (4, 4), (4, 0)]
obstacles = []
sa = SurveyArea(sa_vertices, obstacles)
events = sa.generateEvents()

for e in events:
    print(e)


