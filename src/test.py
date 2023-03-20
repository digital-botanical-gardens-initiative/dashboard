import dash
import plotly.graph_objs as go
import dash_core_components as dcc
from dash import html
from dash.dependencies import Input, Output
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
import base64
from io import BytesIO

from application import app


df = pd.read_csv("/home/mwannier/dashboard/data/230106_frozen_metadata.csv")


graph_component = dcc.Graph(
    id='tsne',
    config={'displayModeBar': False},
    figure={
        'data': [
            go.Scattergl(
                x=df['structure_xlogp'],
                y=df['structure_exact_mass'],
                mode='markers',
                opacity=0.7,
                marker={
                    'size': 5,
                    'color': 'orange',
                    'line': {'width': 0.5, 'color': 'white'}
                },
                name="Decoy"
            )
        ],
        'layout': go.Layout(
            height=400,
            xaxis={'title': 'X'},
            yaxis={'title': 'Y'},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 1, 'y': 1},
            hovermode=False,
            dragmode='select'
        )
    }
)

image_component = html.Img(id="structure-image")

layout = html.Div([
    html.Div([graph_component]),
    html.Div([image_component])
])

@app.callback(
    Output('structure-image', 'src'),
    [Input('tsne', 'selectedData')])
def display_selected_data(selectedData):
    max_structs = 12
    structs_per_row = 6
    empty_plot = "data:image/gif;base64,R0lGODlhAQABAAAAACwAAAAAAQABAAA="
    if selectedData:
        if len(selectedData['points']) == 0:
            return empty_plot
        match_idx = [x['pointIndex'] for x in selectedData['points']]
        smiles_list = [Chem.MolFromSmiles(x) for x in list(df.iloc[match_idx].structure_smiles)]
        name_list = list(df.iloc[match_idx].structure_nameTraditional)
        name_list = [x for x in zip(name_list)]
        img = MolsToGridImage(smiles_list[0:max_structs], molsPerRow=structs_per_row, legends=name_list)
        buffered = BytesIO()
        img.save(buffered, format="JPEG")
        encoded_image = base64.b64encode(buffered.getvalue())
        src_str = 'data:image/png;base64,{}'.format(encoded_image.decode())
    else:
        return empty_plot
    return src_str


dash.register_page(' test', path='/test', layout=layout, icon="bi bi-house")
