import dash
from dash import Dash, html, dcc, Input, Output
import plotly.express as px
import pandas as pd

from application import app

# load a DataFrame with smiles
df_esol = pd.read_csv("/home/mwannier/dashboard/data/230106_frozen_metadata.csv")

piechart = px.pie(df_esol,names='structure_taxonomy_classyfire_02superclass')
piechart.update_traces(textposition='inside')

layout = html.Div(children=[
    html.H3(children='Piechart'),
    dcc.Dropdown(id='names',
        options=['structure_taxonomy_classyfire_02superclass', 'structure_taxonomy_classyfire_03class', 'structure_taxonomy_classyfire_04directparent'],
        clearable=False),
    dcc.Graph(id="graph")
])

@app.callback(
    Output("graph", "figure"), 
    Input("names", "value")) 
def generate_chart(col_name): 
    fig = px.pie(df_esol, names=col_name)
    fig.update_traces(textposition='inside')
    return fig

dash.register_page(' Piecharts', path='/piecharts', layout=layout, icon="bi bi-pie-chart")
