import json
import urllib.request as urlreq
import dash
from dash.dependencies import Input, Output
import dash_bio as dashbio
from dash import html

from application import app

model_data = urlreq.urlopen(
    'https://git.io/mol2d_buckminsterfullerene.json'
).read().decode('utf-8')

model_data = json.loads(model_data)

layout = html.Div([
    dashbio.Molecule2dViewer(
        id='dashbio-default-molecule2d',
        modelData=model_data
    ),
    html.Hr(),
    html.Div(id='default-molecule2d-output')
])

@app.callback(
    Output('default-molecule2d-output', 'children'),
    Input('dashbio-default-molecule2d', 'selectedAtomIds')
)
def update_selected_atoms(ids):
    if ids is None or len(ids) == 0:
        return "No atom has been selected. Select atoms by clicking on them."
    return "Selected atom IDs: {}.".format(', '.join([str(i) for i in ids]))


dash.register_page(' 2D molecules visualization', path='/2D-mol', layout=layout, icon="bi bi-house")
