import dash
from dash import html, dcc, Input, Output, dash_table, State
import dash_bootstrap_components as dbc
from application import app
import pandas as pd
import dash_bio as dashbio
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64
import os
from tqdm import tqdm
from dash.exceptions import PreventUpdate
import time


CSV_PATH = f'{os.getcwd()}/data/230106_frozen_metadata.csv'


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

    # Change smiles to img in html format
    row['structure_smiles'] = smile_to_img_md(row['structure_smiles'])
    row['structure_smiles_2D'] = smile_to_img_md(row['structure_smiles_2D'])

    return row


def apply_md_data(df):
    modified_rows = []
    for index, row in tqdm(df.iterrows(), total=len(df)):
        if row['changed'] == 0:
            modified_row = md_data(row)
            modified_row['changed'] = 1
            modified_rows.append(modified_row)
        else:
            modified_rows.append(row)
    return pd.DataFrame(modified_rows)



header = html.H3('Welcome to home page!')

    # Read the data from the CSV file
df = pd.read_csv(CSV_PATH, nrows=100)
df['changed'] = 0


layout = html.Div(
    [
        header,
        dcc.Input(
            id="input_text",
            type="text",
            placeholder="input type text",
        ),
        html.Button('Submit', id='submit-val', n_clicks=0),
        html.Div(id='progress-container'),
        html.Div(id='datatable-container', style={'height': '100vh',
                                                    'width': '100vw',
                                                    'margin-right': '50',
                                                    'margin-left':'50'})
    ]
)


@app.callback(
    Output('datatable-container', 'children'),
    Output('progress-container', 'children'),
    Input('submit-val', 'n_clicks'),
    State("input_text", "value"),
)
def search_dataframe_vectorized(n_clicks, user_input):

    if user_input is None or user_input == "":
        raise PreventUpdate  # Prevent the callback from running

    # Create the progress bar
    progress_bar = html.Div([
        html.Div('', id='progress-bar', style={'width': '0%'}),
        html.Div('0%', id='progress-label'),
    ])

    # Set the total number of iterations
    total_iterations = 100

    # Create a list to store the progress bar updates
    progress_updates = []

    # Loop through the iterations
    for i in tqdm(range(total_iterations)):
        time.sleep(0.1)  # replace with your function
        progress_updates.append(html.Div('', style={'width': f'{i}%'}))
        progress_label = f'{i+1}%'
        if i == total_iterations-1:
            progress_label = 'Done!'
        progress_updates.append(html.Div(progress_label, style={'margin-left': f'{i}%'}))

    # Use the str.contains method to create a boolean mask of rows that contain the user input
    mask = df.apply(lambda x: x.str.contains(user_input, case=False), axis=1)

    # Filter the DataFrame to keep only the rows that match the mask
    matching_rows = df[mask.any(axis=1)]

    # Create a boolean mask of the rows where 'changed' is equal to 0
    mask_to_change = matching_rows['changed'] == 0

    modified_rows = apply_md_data(matching_rows)

    # Convert the matching rows to a list of dictionaries
    data = modified_rows.to_dict('records')

    # Create a Dash DataTable using the list of dictionaries
    datatable = dash_table.DataTable(
        id='datatable',
        columns=[{'id': col, 'name': col, 'presentation': 'markdown'} if col in ['structure_wikidata','organism_wikidata','reference_wikidata','structure_smiles_2D','structure_smiles', 'reference_doi'] else {'id': col, 'name': col} for col in matching_rows.columns],
        data=data,
        markdown_options={"html": True},
        style_table={'overflowX': 'scroll', 'marginRight': '100'},
        style_cell={
            'whiteSpace': 'normal',
            'textAlign': 'left',
        },
        style_as_list_view=True,
        sort_action='native'
    )

    # Return the DataTable and the progress bar updates
    return datatable, progress_bar

dash.register_page('Home', path='/', layout=layout, icon="bi bi-house")


