import dash
from dash import html, dcc, Input, Output, dash_table, State, callback, register_page, clientside_callback
import pandas as pd
import dash_bio as dashbio
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64
from dash.exceptions import PreventUpdate
import dask.dataframe as dd
from sqlalchemy.orm import sessionmaker
from sqlalchemy import text
from sqlalchemy import create_engine
from dask.delayed import delayed
from rdkit.Chem import AllChem
from rdkit import DataStructs
from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple
import psycopg2
import os
import dash_bootstrap_components as dbc


# import env variable
from dotenv import load_dotenv
load_dotenv()

db_name=os.getenv('DB_NAME')
db_pwd=os.getenv('DB_PWD')
db_user=os.getenv('DB_USR')


META = {'structure_wikidata': 'bool', 'structure_inchikey': 'bool', 'structure_inchi': 'bool', 'structure_smiles': 'bool', 'structure_molecular_formula': 'bool', 'structure_exact_mass': 'float64', 'structure_xlogp': 'float64', 'structure_smiles_2D': 'bool', 'structure_cid': 'float64', 'structure_nameIupac': 'bool', 'structure_nameTraditional': 'bool', 'structure_stereocenters_total': 'float64', 'structure_stereocenters_unspecified': 'float64', 'structure_taxonomy_npclassifier_01pathway': 'bool', 'structure_taxonomy_npclassifier_02superclass': 'bool', 'structure_taxonomy_npclassifier_03class': 'bool', 'structure_taxonomy_classyfire_chemontid': 'float64', 'structure_taxonomy_classyfire_01kingdom': 'bool', 'structure_taxonomy_classyfire_02superclass': 'bool', 'structure_taxonomy_classyfire_03class': 'bool', 'structure_taxonomy_classyfire_04directparent': 'bool', 'organism_wikidata': 'bool', 'organism_name': 'bool', 'organism_taxonomy_gbifid': 'bool', 'organism_taxonomy_ncbiid': 'float64', 'organism_taxonomy_ottid': 'float64', 'organism_taxonomy_01domain': 'bool', 'organism_taxonomy_02kingdom': 'bool', 'organism_taxonomy_03phylum': 'bool', 'organism_taxonomy_04class': 'bool', 'organism_taxonomy_05order': 'bool', 'organism_taxonomy_06family': 'bool', 'organism_taxonomy_07tribe': 'bool', 'organism_taxonomy_08genus': 'bool', 'organism_taxonomy_09species': 'bool', 'organism_taxonomy_10varietas': 'bool', 'reference_wikidata': 'bool', 'reference_doi': 'bool', 'manual_validation': 'bool', 'changed': 'float64'}

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
    row['structure_nameTraditional'] = '['+ row['structure_nameTraditional'] + '](/element/' + row['structure_nameTraditional'] + ')'
    row['organism_name'] = '['+ row['organism_name'] + '](/organism/' + row['organism_name'] + ')'


    # Change smiles to img in html format
    row['structure_smiles'] = smile_to_img_md(row['structure_smiles'])
    row['structure_smiles_2D'] = smile_to_img_md(row['structure_smiles_2D'])

    return row


def apply_md_data(ddf):
    modified_ddf = ddf.copy()
    modified_ddf = modified_ddf.apply(md_data, axis=1, meta=META)
    return modified_ddf

def get_molecules_from_db(batch_size: int = 1000) -> List[Tuple[int, str]]:
    connection = psycopg2.connect(f"dbname={db_name} user={db_user} password={db_pwd} host=localhost")
    cursor = connection.cursor()

    cursor.execute(f"SELECT COUNT(*) FROM {db_name}")
    total_molecules = cursor.fetchone()[0]

    for offset in range(0, total_molecules, batch_size):
        cursor.execute(f"SELECT structure_nameTraditional, structure_smiles_2D FROM {db_name} LIMIT {batch_size} OFFSET {offset}")
        yield cursor.fetchall()

    cursor.close()
    connection.close()

def calculate_similarity(reference_fp, mol_smiles: str, threshold: float) -> Tuple[str, float]:
    mol = Chem.MolFromSmiles(mol_smiles)
    mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    tanimoto_similarity = DataStructs.TanimotoSimilarity(reference_fp, mol_fp)

    if tanimoto_similarity >= threshold:
        return mol_smiles, tanimoto_similarity
    return None

def find_similar_molecules(reference_smiles: str, threshold: float, n_workers: int = 4) -> List[Tuple[str, float]]:
    reference_mol = Chem.MolFromSmiles(reference_smiles)
    reference_fp = AllChem.GetMorganFingerprintAsBitVect(reference_mol, 2, nBits=2048)

    similar_molecules = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        for batch in get_molecules_from_db():
            future_results = [executor.submit(calculate_similarity, reference_fp, mol_smiles, threshold) for _, mol_smiles in batch]
            for future in future_results:
                result = future.result()
                if result:
                    similar_molecules.append(result[0])  # Append only the SMILES string

    return similar_molecules

register_page(__name__, path="/explore")

header = html.H3('Explore data:')

# Replace the following with your own database connection information
db_connection_str = f'postgresql://{db_user}:{db_pwd}@localhost/{db_name}'
engine = create_engine(db_connection_str)
Session = sessionmaker(bind=engine)

# Replace 'your_table_name' with the actual table name in your database
query = f"SELECT * FROM {db_name}"

# Define a function to read chunks of data from the database
@delayed
def load_chunk(chunk_id, chunk_size):
    offset = chunk_id * chunk_size
    chunk_query = f"{query} LIMIT {chunk_size} OFFSET {offset}"
    with engine.connect() as conn:
        return pd.read_sql(text(chunk_query), con=conn)

# Set the number of chunks and chunk size based on your dataset size and available memory
num_chunks = 10
chunk_size = 10000

# Read data from the PostgreSQL table into Dask DataFrame using the delayed function
chunks = [load_chunk(chunk_id, chunk_size) for chunk_id in range(num_chunks)]
ddf = dd.from_delayed(chunks)

ddf['changed'] = 0


layout = html.Div(
    [   header,
        html.Hr(),
        dcc.Dropdown(options=[{'label': 'By name', 'value': 'by_name'}, {'label': 'By structure', 'value': 'by_structure'}], id='dropdown'),
        html.Br(),
        dcc.Loading(id='loading', children=[
                                            html.Div(id='search_layout'),
                                            html.Div(id='hidden-data-text', style={'display': 'none'}),
                                            html.Div(id='hidden-data-structure', style={'display': 'none'}),
                                            html.Div(id='hidden-smiles'),  # Add this line
                                            html.Div(id='datatable-container', style={'height': '100vh',
                                                                                        'width': '100vw',
                                                                                        'margin-right': '50',
                                                                                        'margin-left':'50'})],
                    type ="circle"
        )
])


@callback(
    Output('search_layout', 'children'),
    Input('dropdown', 'value'),
    prevent_initial_callbacks=True
)
def search_layout(search_type):
    if search_type == 'by_name':
        layout = html.Div([
                        dcc.Input(
                            id="input_text",
                            type="text",
                            placeholder="input type text",
                        ),
                        dcc.Checklist(options=[{'label': 'Exact', 'value': 'exact'}], value=[], id='exact-search-name'),
                        html.Button(
                            'Submit',
                            id='submit-val', 
                            n_clicks=0)])
    elif search_type == 'by_structure':
        layout = html.Div([html.Div(id='molecule-editor', style={'width': '800px', 'height': '600px'}),
                            html.Br(),
                            dbc.Input(id='smiles-input', type='text', placeholder='SMILES value will be displayed here...')])
    else:
        raise PreventUpdate

    return layout

@callback(
    Output('smiles-output', 'children'),
    Input('smiles-input', 'value')
)
def update_smiles_output(smiles):
    return smiles



# Create a callback for the 'by text' case
@callback(
    Output('hidden-data-text', 'children'),
    Input('submit-val', 'n_clicks'),
    State("input_text", "value"),
    prevent_initial_callbacks=True
)
def search_text(n_clicks_text, user_input):
    if n_clicks_text is None or user_input is None or user_input == "":
        return None
    if user_input is not None and user_input != "" and n_clicks_text is not None and n_clicks_text > 0:
        # Use the str.contains method to create a boolean mask of rows that contain the user input
        mask_text = ddf.apply(lambda x: x.str.contains(user_input, case=False), axis=1, meta=META)

        # Filter the DataFrame to keep only the rows that match the mask
        matching_rows_text = ddf[mask_text.any(axis=1)]

        modified_rows_text = apply_md_data(matching_rows_text)

        return modified_rows_text.compute().to_json(date_format='iso', orient='split')
    
# Create a callback for the 'by structure' case
@callback(
    Output('hidden-data-structure', 'children'),
    Input('submit-structure', 'n_clicks'),
    Input('exact-search-structure', 'value'),
    Input("input_structure", "smiles"),
    prevent_initial_callbacks=True
)
def search_structure(n_clicks_structure, exact_search_value, value_structure):
    print(value_structure)

    if n_clicks_structure is None or value_structure is None:
        return None
    
    if value_structure is not None and n_clicks_structure is not None and n_clicks_structure > 0:
        smiles = value_structure

        if exact_search_value:
            mask_structure = ddf['structure_smiles_2D'] == smiles
            matching_rows_structure = ddf[mask_structure]
        else:
            similar_molecules = find_similar_molecules(smiles, 0.85)
            similar_smiles = [smiles for smiles, _ in similar_molecules]
            matching_rows_structure = ddf[ddf['structure_smiles_2D'].isin(similar_smiles)]


        modified_rows_structure = apply_md_data(matching_rows_structure)

    return modified_rows_structure.compute().to_json(date_format='iso', orient='split')


@callback(
    Output('datatable-container', 'children'),
    Input('hidden-data-text', 'children'),
    Input('hidden-data-structure', 'children'),
    Input('dropdown', 'value'),
    prevent_initial_callbacks=True
)
def display_datatable(data_text_json, data_structure_json, dropdown_value):
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

    # If the dropdown is the trigger, clear the datatable
    if triggered_id == 'dropdown':
        return []

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
        columns=[{'id': col, 'name': col, 'presentation': 'markdown'} if col in ['organism_name', 'structure_wikidata', 'organism_wikidata', 'reference_wikidata', 'structure_smiles_2D', 'structure_smiles', 'reference_doi', 'structure_nameTraditional'] else {'id': col, 'name': col} for col in ddf.columns],
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
        virtualization=True,  # Enable virtualization for lazy loading
        page_size=10  # Set the number of rows per page
    )





