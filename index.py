import dash
from dash import dcc, html, Output, Input, State
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.exceptions import PreventUpdate

from src import mol_visu, home, smiles_visu, molplotly_layout, piechart_layout, element_visu, organisms#, tmap_visu
#from GNPS_LCMSDashboard import layout

from application import app

home_page = ['Home']
    
sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div([
                    html.I(className='bi bi-house'),
                    ' Home'
                ], className="ms-2"),
            ],
            href='/',
            active="exact",
        ),
        dmc.Accordion(
        dmc.AccordionItem(
            [
            dbc.NavLink(
                [
                    html.Div([
                        html.I(className=page["icon"]),
                        page["name"]
                    ], className="ms-2"),
                ],
                href=page["path"],
                active="exact",
            )
            for page in dash.page_registry.values() if page['name'] not in home_page

            ], label='Explore')
        )
    ],
    vertical=True,
    pills=True,
    className="bg-light"
)



app.layout = dbc.Container([
    dcc.Location(id="url"),
    dbc.Offcanvas([
        sidebar,
    ],
        id="offcanvas",
        title="Welcome to .../",
        is_open=False,
    ),
    html.Div([
        dbc.Button(html.I(className="bi bi-list-ul"),  # "Menu",
                   id="open-offcanvas", n_clicks=0),
    ]),

    dbc.Row([
        dbc.Col(html.Div("Digital Botanical Garden Initiative",
                         style={'fontSize': 50, 'textAlign': 'center', 'font-family':'Adorable'}))
    ]),

    html.Hr(),

    dbc.Row([
        dbc.Col([
            dash.page_container
        ], xs=8, sm=8, md=10, lg=10, xl=10, xxl=10)
    ]),
    html.Div(id='page-content')
     ], fluid=True
#style={'background-image':'url("/assets/cactus.jpg")', 'background-repeat': 'no-repeat','opacity':'0.4'}
 )


@app.callback(Output('page-content', 'children'),
              Input('url', 'pathname'))
def display_page(pathname):
    if pathname.startswith('/element/'):
        element_name = pathname.split('/')[-1]
        return element_visu.element_visualization(element_name, 'structure_nameTraditional')
    elif pathname.startswith('/organism/'):
        organism_name = pathname.split('/')[-1]
        return organisms.organisms_visu(organism_name)


@app.callback(
    Output('offcanvas', 'is_open'),
    Input('open-offcanvas','n_clicks'),
    State('offcanvas', 'is_open'),
)
def openNav(n1, o):
    if n1:
        return not o

if __name__ == "__main__":
    app.run_server(debug=True)