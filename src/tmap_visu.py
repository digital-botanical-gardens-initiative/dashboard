import tmap as tm
import numpy as np
from matplotlib import pyplot as plt
import dash
from dash import html
import chart_studio.plotly as py

import plotly.express as px



n = 25
edge_list = []
# Create a random graph
for i in range(n):
    for j in np.random.randint(0, high=n, size=2):
        edge_list.append([i, j, np.random.rand(1)])
# Compute the layout
x, y, s, t, _ = tm.layout_from_edge_list(
    n, edge_list, create_mst=False
)
# Plot the edges
for i in range(len(s)):
    plt.plot(
        [x[s[i]], x[t[i]]],
        [y[s[i]], y[t[i]]],
        "k-",
        linewidth=0.5,
        alpha=0.5,
        zorder=1,
    )
# Plot the vertices
plt.scatter(x, y, zorder=2)
plt.tight_layout()
plt.show()

tmap_plot = py.plot_mpl(plt)

layout = html.Div(children=[
    tmap_plot
])

#dash.register_page(' Tmap', path='/tmap', layout=layout, icon="bi bi-arrow-through-heart-fill")
