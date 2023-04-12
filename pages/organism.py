import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from dash import html, dcc, Input, Output, register_page, callback
from io import BytesIO
import psycopg2
import pandas as pd
import os


# import env variable
from dotenv import load_dotenv
load_dotenv()

db_name=os.getenv('DB_NAME')
db_pwd=os.getenv('DB_PWD')
db_user=os.getenv('DB_USR')


register_page(
    __name__,
    path="/organism"
)


connection = psycopg2.connect(
        dbname=db_name,
        user=db_user,
        password=db_pwd,
        host="localhost",
    )
query = f'SELECT organism_name FROM {db_name}'
df = pd.read_sql_query(query, connection)

# Close the connection
connection.close()
organisms = sorted(list(set(df['organism_name'].astype(str))))

# Generate alphabet buttons
alphabet = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
alphabet_buttons = dbc.ButtonGroup([dbc.Button(letter, id=letter, color='primary') for letter in alphabet])


layout = dbc.Container([
    html.H1("Species List"),
    dbc.Row([
        dbc.Col(alphabet_buttons)
    ]),
    html.Div(id='species-list')
], fluid=True)

@callback(
    Output('species-list', 'children'),
    [Input(letter, 'n_clicks') for letter in alphabet]
)
def update_species_list(*args):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    filtered_species = [s for s in organisms if s[0].upper() == button_id]

    num_columns = 3
    column_length = len(filtered_species) // num_columns + (len(filtered_species) % num_columns > 0)
    columns = [filtered_species[i:i + column_length] for i in range(0, len(filtered_species), column_length)]

    return dbc.Row([
                    dbc.Col([
                            html.Ul([html.Li(dcc.Link(item,href=f'/organism/{item}')) for item in col])
                ]) for col in columns])

