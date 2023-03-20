import dash
from dash import html, dcc, Input, Output, dash_table, State
from application import app
import pandas as pd
import dash_bio as dashbio
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import base64

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

def md_data(rows_masked):
    unchanged_mask = (rows_masked['changed'] == 0)
    
    #Display data as link in md format
    rows_masked.loc[unchanged_mask,'structure_wikidata'] = rows_masked.loc[unchanged_mask,'structure_wikidata'].apply(lambda x: '['+ x.split('/')[-1] + '](' + x + ')')
    rows_masked.loc[unchanged_mask,'organism_wikidata'] = rows_masked.loc[unchanged_mask,'organism_wikidata'].apply(lambda x: '['+ x.split('/')[-1] + '](' + x + ')')
    rows_masked.loc[unchanged_mask,'reference_wikidata'] = rows_masked.loc[unchanged_mask,'reference_wikidata'].apply(lambda x: '['+ x.split('/')[-1] + '](' + x + ')')
    rows_masked.loc[unchanged_mask, 'reference_doi'] = rows_masked.loc[unchanged_mask, 'reference_doi'].apply(lambda x: '[' + x + '](https://doi.org/' + x + ')')

    # Change smiles to img in html format
    rows_masked.loc[unchanged_mask, 'structure_smiles'] = rows_masked.loc[unchanged_mask,'structure_smiles'].apply(lambda x: smile_to_img_md(x))
    rows_masked.loc[unchanged_mask,'structure_smiles_2D'] = rows_masked.loc[unchanged_mask,'structure_smiles_2D'].apply(lambda x: smile_to_img_md(x))

    #df.loc[unchanged_mask, 'changed'] = 1




header = html.H3('Welcome to home page!')

    # Read the data from the CSV file
df = pd.read_csv("/home/mwannier/dashboard/data/230106_frozen_metadata.csv")
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
        html.Div(id='datatable-container', style={'height': '100vh',
                                                    'width': '100vw',
                                                    'margin-right': '50',
                                                    'margin-left':'50'})
    ]
)


@app.callback(
    Output('datatable-container', 'children'),
    Input('submit-val', 'n_clicks'),
    State("input_text", "value"),
)
def search_dataframe_vectorized(n_clicks, user_input):

    if user_input is None or user_input == "":
        return html.Div()
    # Otherwise, return the DataTable
    else:
        # Use the str.contains method to create a boolean mask of rows that contain the user input
        mask = df.apply(lambda x: x.str.contains(user_input, case=False), axis=1)

        # Filter the DataFrame to keep only the rows that match the mask
        matching_rows = df[mask.any(axis=1)]

        md_data(matching_rows)

        # Convert the matching rows to a list of dictionaries
        data = matching_rows.to_dict('records')

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

            
        return datatable


dash.register_page('Home', path='/', layout=layout, icon="bi bi-house")


