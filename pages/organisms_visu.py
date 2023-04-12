from urllib.parse import unquote
from dash import html, dcc, dash_table, register_page
import dash_bootstrap_components as dbc
from urllib.parse import unquote
import psycopg2
import pandas as pd
import os

# import env variable
from dotenv import load_dotenv
load_dotenv()

db_name=os.getenv('DB_NAME')
db_pwd=os.getenv('DB_PWD')
db_user=os.getenv('DB_USR')


# register page
register_page(
    __name__,
    path_template="/organism/<organismname>",
    path="/organism/Acarnus erithacus",
)


def layout(organismname=None, **other_unknown_query_strings):
    organismname = unquote(organismname)
    connection = psycopg2.connect(
        dbname=db_name,
        user=db_user,
        password=db_pwd,
        host="localhost",
    )
    query = f'SELECT * FROM {db_name} WHERE "organism_name" = \'{organismname}\''
    df = pd.read_sql_query(query, connection)
    # Close the connection
    connection.close()

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

    elements = sorted(list(set(df['structure_nameTraditional'].astype(str))))
    num_columns = 2
    column_length = len(elements) // num_columns + (len(elements) % num_columns > 0)
    columns = [elements[i:i + column_length] for i in range(0, len(elements), column_length)]

    layout = html.Div([                html.Div([                        html.H2([organismname, ' (', html.A(wikidata_id, href=wikidata),')'], style={'font-family': 'fantasy', 'text-align': 'center'}),
                        html.Br(),
                        html.Hr(),
                        html.Div([html.H4('Taxonomy'), datatable]),
                        html.Hr(),
                        html.Div([html.H4('Molecules'),                                    dbc.Row([                                        dbc.Col([                                            html.Ul([html.Li(dcc.Link(item,href=f'/element/{item}')) for item in col])
                                        ]) for col in columns
                                    ])
                            ])
                ])
    ], style={'height': '100vh',
                'width': '100vw',
                'margin-right': '50',
                'margin-left':'50'})

    return layout