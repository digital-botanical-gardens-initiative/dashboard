import dash
from dash import html, dcc, Input, Output, dash_table, State
import dash_bootstrap_components as dbc
from application import app
import dash_bio as dashbio
import dask.dataframe as dd
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64
import os



def smile_to_img(smile):
    buffered = BytesIO()
    d2d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    opts = d2d.drawOptions()
    opts.clearBackground = False
    for i in smile:
        d2d.DrawMolecule(Chem.MolFromSmiles(i))
        d2d.FinishDrawing()
        img_str = d2d.GetDrawingText()
        buffered.write(str.encode(img_str))
        img_str = base64.b64encode(buffered.getvalue())
        img_str = f"data:image/svg+xml;base64,{repr(img_str)[2:-1]}"
    return img_str

def element_visualization(element,column):
    CSV_PATH = f'{os.getcwd()}/data/230106_frozen_metadata.csv'

    df = dd.read_csv(CSV_PATH, dtype={'manual_validation': 'object',
                                        'organism_taxonomy_07tribe': 'object',
                                        'organism_taxonomy_10varietas': 'object',
                                        'organism_taxonomy_gbifid': 'object',
                                        'organism_taxonomy_ncbiid': 'float64',
                                        'organism_taxonomy_ottid': 'float64',
                                        'structure_cid': 'float64',
                                        'structure_taxonomy_classyfire_chemontid': 'float64'})
    df = df[(df[column] == element)]
    df = df.compute()

    name = list(set(df['structure_nameTraditional'].astype(str)))[0]
    formula = list(set(df['structure_molecular_formula'].astype(str)))[0]
    wikidata = list(set(df['structure_wikidata'].astype(str)))[0]
    wikidata_id = wikidata.split('/')[-1]
    smile = list(set(df['structure_smiles'].astype(str)))
    smile_2d = list(set(df['structure_smiles_2D'].astype(str)))

    img_2d_str = smile_to_img(smile_2d)
    img_str = smile_to_img(smile)

    phylogeny = df.iloc[0, 13:20].copy().to_frame().reset_index()
    datatable = dash_table.DataTable(
        id='datatable',
        columns=[{"name": str(i), "id": str(i)} for i in phylogeny.columns],
        data=phylogeny.to_dict('records'),
        markdown_options={"html": True},
        style_table={'overflowX': 'scroll', 'marginRight': '100'},
        style_cell={
            'whiteSpace': 'normal',
            'textAlign': 'left',
        },
        style_as_list_view=True,
        sort_action='native'
    )

    organisms = sorted(list(set(df['organism_name'].astype(str))))
    #organisms_md = '\n'.join([f'- {item}' for item in organisms])
    num_columns = 4
    column_length = len(organisms) // num_columns + (len(organisms) % num_columns > 0)
    columns = [organisms[i:i + column_length] for i in range(0, len(organisms), column_length)]


    layout = html.Div([
                html.Div([
                        html.H2([name, ' (', html.A(wikidata_id, href=wikidata),')'], style={'font-family': 'fantasy', 'text-align': 'center'}),
                        html.Br(),
                        dbc.Row([
                                html.H4('Chemical structures'),
                                html.H6(f'molecular formula : {formula}'),
                                dbc.Col([
                                html.Img(
                                src=img_str, 
                                title='smile', 
                                style={
                                    "width": "50%",
                                    "background-color": f"rgba(255,255,255,{0.7})",
                                    'display': 'block',
                                    'margin-left': 'auto',
                                    'margin-right': 'auto'
                                }),
                                html.Caption('3D')]),
                                dbc.Col([
                                html.Img(src=img_2d_str,
                                title='smile 2D',      
                                style={
                                    "width": "50%",
                                    "background-color": f"rgba(255,255,255,{0.7})",
                                    'display': 'block',
                                    'margin-left': 'auto',
                                    'margin-right': 'auto'
                                }),
                                html.Caption('2D')])
                        ]),
                        html.Hr(),
                        html.Div([html.H4('Phylogeny'), datatable]),
                        html.Hr(),
                        html.Div([html.H4('Organisms'),
                                    dbc.Row([
                                        dbc.Col([
                                            html.Ul([html.Li(dcc.Link(item,href=f'/organism/{item}')) for item in col])
                                        ]) for col in columns
                                    ])
                            ])
                ])
    ], style={'height': '100vh',
                'width': '100vw',
                'margin-right': '50',
                'margin-left':'50'})

    return layout


#layout = element_visualization('http://www.wikidata.org/entity/Q43656','structure_wikidata')

#dash.register_page('element', path='/element', layout=layout, icon="bi bi-house")
