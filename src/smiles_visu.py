from dash import Dash, html
import dash_bio as dashbio
import dash


layout = html.Div([
    dashbio.Jsme(
        smiles='O=C(Nc1cccc(Cl)c1)c3cncc4nnc(c2ccc(OC(F)F)cc2)n34',
    ),
])

dash.register_page('Smiles visualization', path='/smiles', layout=layout, icon="bi bi-house")
