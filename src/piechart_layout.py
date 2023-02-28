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
    html.H1(children='DBGI Dashboard'),
    html.Div(children='''Dash: A web application framework for your data.'''),
    dcc.Graph(id="graph"),
    dcc.Dropdown(id='names',
        options=['structure_taxonomy_classyfire_02superclass', 'structure_taxonomy_classyfire_03class', 'structure_taxonomy_classyfire_04directparent'],
        clearable=False)
])

@app.callback(
    Output("graph", "children"), 
    Input("names", "value")) 
def generate_chart(value): 
    fig = px.pie(df_esol, names=value)
    fig.update_traces(textposition='inside')
    return fig

dash.register_page(' Page 2', path='/page-2', layout=layout, icon="bi bi-pie-chart")
