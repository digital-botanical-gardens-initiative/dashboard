import dash
from dash import html, dcc, Input, Output, dash_table, State, callback, register_page
import dash_bootstrap_components as dbc
import pandas as pd
import dash_bio as dashbio
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64
import os
from dash.exceptions import PreventUpdate

def smile_to_img_md(smile):
    buffered = BytesIO()
    d2d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    opts = d2d.drawOptions()
    opts.clearBackground = False
    d2d.DrawMolecule(Chem.MolFromSmiles(smile))
    d2d.FinishDrawing()
    img_str = d2d.GetDrawingText()
    buffered.write(str.encode(img_str))
    img_str = base64.b64encode(buffered.getvalue())
    img_str = f"data:image/svg+xml;base64,{repr(img_str)[2:-1]}"
    img = "<img src='" + img_str + "' height='90' />"
    return img

def md_data(row):
    #Display data as link in md format
    row['structure_wikidata'] = '['+ row['structure_wikidata'].split('/')[-1] + '](' + row['structure_wikidata'] + ')'
    row['organism_wikidata'] = '['+ row['organism_wikidata'].split('/')[-1] + '](' + row['organism_wikidata'] + ')'
    row['reference_wikidata'] = '['+ row['reference_wikidata'].split('/')[-1] + '](' + row['reference_wikidata'] + ')'
    row['reference_doi'] = '[' + row['reference_doi'] + '](https://doi.org/' + row['reference_doi'] + ')'
    row['structure_nameTraditional'] = '['+ row['structure_nameTraditional'].split('/')[-1] + '](/element/' + row['structure_nameTraditional'] + ')'
    row['organism_name'] = '['+ row['organism_name'].split('/')[-1] + '](/organism/' + row['organism_name'] + ')'
    

    # Change smiles to img in html format
    row['structure_smiles'] = smile_to_img_md(row['structure_smiles'])
    row['structure_smiles_2D'] = smile_to_img_md(row['structure_smiles_2D'])

    return row


def apply_md_data(df):
    modified_df = df.copy()
    modified_df = modified_df.apply(md_data, axis=1)
    return modified_df

register_page(__name__, path="/explore")

CSV_PATH = f'{os.getcwd()}/data/230106_frozen_metadata.csv'

header = html.H3('Explore data:')
    # Read the data from the CSV file
df = pd.read_csv(CSV_PATH, dtype={'manual_validation': 'object',
                                        'organism_taxonomy_07tribe': 'object',
                                        'organism_taxonomy_10varietas': 'object',
                                        'organism_taxonomy_gbifid': 'object',
                                        'organism_taxonomy_ncbiid': 'float64',
                                        'organism_taxonomy_ottid': 'float64',
                                        'structure_cid': 'float64',
                                        'structure_taxonomy_classyfire_chemontid': 'float64'},nrows=10000)
df['changed'] = 0


layout = html.Div(
    [   header,
        html.Hr(),
        dcc.Dropdown(['by name','by structure'],id='dropdown'),
        html.Br(),
        dcc.Loading(id='loading', children=[
                                            html.Div(id='search_layout'),
                                            html.Div(id='hidden-data-text', style={'display': 'none'}),
                                            html.Div(id='hidden-data-structure', style={'display': 'none'}),
                                            html.Div(id='datatable-container', style={'height': '100vh',
                                                                                        'width': '100vw',
                                                                                        'margin-right': '50',
                                                                                        'margin-left':'50'})],
                    type ="circle"
        )
])


@callback(
    Output('search_layout', 'children'),
    Input('dropdown', 'value')
)
def search_layout(search_type):
    if search_type == 'by name':
        layout = html.Div([
                        dcc.Input(
                            id="input_text",
                            type="text",
                            placeholder="input type text",
                        ),
                        dcc.Checklist(['Exact'],[], id='exact-search'),
                        html.Button(
                            'Submit',
                            id='submit-val', 
                            n_clicks=0)])
    elif search_type == 'by structure':
        layout = html.Div([
            dashbio.Jsme( 
                id='input_structure',
                width='70%',
                height='50vh',
            ),
            dcc.Checklist(['Exact'],[], id='exact-search'),
            html.Button(
                'Submit',
                id='submit-structure',
                n_clicks=0
            )
        ])
    else:
        raise PreventUpdate

    return layout

# Create a callback for the 'by text' case
@callback(
    Output('hidden-data-text', 'children'),
    Input('submit-val', 'n_clicks'),
    State("input_text", "value")
)
def search_text(n_clicks_text, user_input):
    if n_clicks_text is None or user_input is None or user_input == "":
        return None
    if user_input is not None and user_input != "" and n_clicks_text is not None and n_clicks_text > 0:
        # Use the str.contains method to create a boolean mask of rows that contain the user input
        mask_text = df.apply(lambda x: x.str.contains(user_input, case=False), axis=1)

        # Filter the DataFrame to keep only the rows that match the mask
        matching_rows_text = df[mask_text.any(axis=1)]

        modified_rows_text = apply_md_data(matching_rows_text)

        return modified_rows_text.to_json(date_format='iso', orient='split')
    
# Create a callback for the 'by structure' case
@callback(
    Output('hidden-data-structure', 'children'),
    Input('submit-structure', 'n_clicks'),
    State("input_structure", "value")
)
def search_structure(n_clicks_structure, value_structure):
    if n_clicks_structure is None or value_structure is None:
        return None
    
    if value_structure is not None and n_clicks_structure is not None and n_clicks_structure > 0:
        smiles = value_structure.getSmiles()
        mask_structure = df['structure_smiles'].apply(lambda x: Chem.MolFromSmiles(x) is not None and Chem.MolFromSmiles(x).GetMol() is not None and Chem.MolFromSmiles(x).GetMol().HasSubstructMatch(Chem.MolFromSmiles(smiles)))
        matching_rows_structure = df[mask_structure]
    
        modified_rows_structure = apply_md_data(matching_rows_structure)

    return modified_rows_structure.to_json(date_format='iso', orient='split')

# Create a callback for displaying the datatable
@callback(
    Output('datatable-container', 'children'),
    Input('hidden-data-text', 'children'),
    Input('hidden-data-structure', 'children'),
)
def display_datatable(data_text_json, data_structure_json):
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if triggered_id == 'hidden-data-text':
        if data_text_json is None:
            raise PreventUpdate
        data_text = pd.read_json(data_text_json, orient='split').to_dict('records')
        datatable_text = create_datatable(data_text)
        return datatable_text
    elif triggered_id == 'hidden-data-structure':
        if data_structure_json is None:
            raise PreventUpdate
        data_structure = pd.read_json(data_structure_json, orient='split').to_dict('records')
        datatable_structure = create_datatable(data_structure)
        return datatable_structure
    else:
        raise PreventUpdate

def create_datatable(data):
    return dash_table.DataTable(
        id='datatable',
        columns=[{'id': col, 'name': col, 'presentation': 'markdown'} if col in ['organism_name', 'structure_wikidata', 'organism_wikidata', 'reference_wikidata', 'structure_smiles_2D', 'structure_smiles', 'reference_doi', 'structure_nameTraditional'] else {'id': col, 'name': col} for col in df.columns],
        data=data,
        markdown_options={"html": True},
        style_table={'overflowX': 'scroll', 'marginRight': '100'},
        style_cell={
            'whiteSpace': 'normal',
            'textAlign': 'left',
        },
        style_as_list_view=True,
        sort_action='native',
        page_action='native',  # Enable pagination
        page_size=10  # Set the number of rows per page
    )



