from dash import Dash, html, dcc, Input, Output
import plotly.express as px
import pandas as pd
from app import app
import vaex


# load a DataFrame with smiles
df_esol = pd.read_csv("/home/mwannier/dbgi_dashboard/data/230106_frozen_metadata.csv")

piechart = px.pie(df_esol,names='structure_taxonomy_classyfire_02superclass')
piechart.update_traces(textposition='inside')

layout = html.Div(children=[
    html.H1(children='DBGI Dashboard'),

    html.Div(children='''
        Dash: A web application framework for your data.
    '''),
    dcc.Graph(id="graph"),
    dcc.Dropdown(id='names',
        options=['structure_taxonomy_classyfire_02superclass', 'structure_taxonomy_classyfire_03class', 'structure_taxonomy_classyfire_04directparent'],
        clearable=False),
])

@app.callback(
    Output("graph", "figure"), 
    Input("names", "value")) 
def generate_chart(names): 
    fig = px.pie(df_esol, names=names)
    fig.update_traces(textposition='inside')
    return fig
