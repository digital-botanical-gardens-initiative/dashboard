from dash import html, dcc, Output, Input
import dash_bio as dashbio
import dash
import pandas as pd

from application import app

df = pd.read_csv("/home/mwannier/dashboard/data/230106_frozen_metadata.csv", nrows=1000)

smiles_unique = df[['structure_taxonomy_classyfire_04directparent','structure_smiles_2D']].drop_duplicates()


layout = html.Div([
    dcc.Dropdown(id='names', options=[{'label': row['structure_taxonomy_classyfire_04directparent'], 
                                       'value': row['structure_smiles_2D']} for index, row in df.iterrows()]),
    html.Div(id='jsme-container')
])

# Define the callback function
@app.callback(Output('jsme-container', 'children'), [Input('names', 'value')])
def update_jsme(value):
    if value:
        return dashbio.Jsme(smiles=value)
    else:
        return ''

dash.register_page(' Smiles visualization', path='/smiles', layout=layout, icon="bi bi-house")
