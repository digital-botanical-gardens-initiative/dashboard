from urllib.parse import unquote
import dash
from dash import html, dcc, Input, Output, dash_table, State, register_page
import dash_bootstrap_components as dbc
import dash_bio as dashbio
import dask.dataframe as dd
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64
import os
from urllib.parse import unquote

register_page(
    __name__,
    path_template="/organism/<organismname>",
    path="/organism/Acarnus erithacus",
)


def layout(organismname=None, **other_unknown_query_strings):
    organismname = unquote(organismname)
    CSV_PATH = f'{os.getcwd()}/data/230106_frozen_metadata.csv'

    df = dd.read_csv(CSV_PATH, dtype={'manual_validation': 'object',
                                        'organism_taxonomy_07tribe': 'object',
                                        'organism_taxonomy_10varietas': 'object',
                                        'organism_taxonomy_gbifid': 'object',
                                        'organism_taxonomy_ncbiid': 'float64',
                                        'organism_taxonomy_ottid': 'float64',
                                        'structure_cid': 'float64',
                                        'structure_taxonomy_classyfire_chemontid': 'float64'})
    
    df = df[(df['organism_name'] == organismname)]
    df = df.compute()

    taxonomy = df.iloc[0, 26:35].copy().to_frame().reset_index()
    datatable = dash_table.DataTable(
        id='datatable',
        columns=[{"name": str(i), "id": str(i)} for i in taxonomy.columns],
        data=taxonomy.to_dict('records'),
        markdown_options={"html": True},
        style_table={'overflowX': 'scroll', 'marginRight': '100'},
        style_cell={
            'whiteSpace': 'normal',
            'textAlign': 'left',
        },
        style_as_list_view=True,
        sort_action='native'
    )


    wikidata = list(set(df['organism_wikidata'].astype(str)))[0]
    wikidata_id = wikidata.split('/')[-1]

    organisms = sorted(list(set(df['structure_nameTraditional'].astype(str))))
    num_columns = 2
    column_length = len(organisms) // num_columns + (len(organisms) % num_columns > 0)
    columns = [organisms[i:i + column_length] for i in range(0, len(organisms), column_length)]

    layout = html.Div([
                html.Div([
                        html.H2([organismname, ' (', html.A(wikidata_id, href=wikidata),')'], style={'font-family': 'fantasy', 'text-align': 'center'}),
                        html.Br(),
                        html.Hr(),
                        html.Div([html.H4('Taxonomy'), datatable]),
                        html.Hr(),
                        html.Div([html.H4('Molecules'),
                                    dbc.Row([
                                        dbc.Col([
                                            html.Ul([html.Li(dcc.Link(item,href=f'/element/{item}')) for item in col])
                                        ]) for col in columns
                                    ])
                            ])
                ])
    ], style={'height': '100vh',
                'width': '100vw',
                'margin-right': '50',
                'margin-left':'50'})

    return layout