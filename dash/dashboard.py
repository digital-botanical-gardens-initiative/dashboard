import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
import webbrowser


app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.config.suppress_callback_exceptions = True


navbar = dbc.NavbarSimple(children=[
        dbc.NavItem(dbc.NavLink("Home", href="/")),
        dbc.DropdownMenu([
            dbc.DropdownMenuItem('Whole dataset',href='/explore'),
            dbc.DropdownMenuItem('Molecules', href='/element'),
            dbc.DropdownMenuItem('Organisms', href='/organism')
        ],
        nav=True,
        in_navbar=True,
        label='Explore'),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("About DBGI", href="/about"),
                dbc.DropdownMenuItem("Documentation", href="/documentation"),
                dbc.DropdownMenuItem('Github',href='https://github.com/digital-botanical-gardens-initiative', target="_blank")
            ],
            nav=True,
            in_navbar=True,
            label="Information",
        ),
    ],
    brand="Digital Botanical Garden Initiative",
    color="green",
    dark=True,
    className="mb-2",
)
'''
navbar = dbc.Navbar(
    dbc.Container([
        html.A(
            dbc.Row([
                dbc.Col(html.Img(src='/assets/logo.png', height='30px')),
                dbc.Col(dbc.NavbarBrand('Digital Botanical Garden Initiative', className='ms-2')),
                ],
                align='right',
                className='g-0',
            ),
            href='/',
            style={"textDecoration":"none"}
        ),
        dbc.NavItem(dbc.NavLink("Home", href="/")),
        dbc.NavItem(dbc.NavLink('Explore', href='/explore')),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("About DBGI", href="/about"),
                dbc.DropdownMenuItem("Documentation", href="/documentation"),
                dbc.DropdownMenuItem('Github',href='https://github.com/digital-botanical-gardens-initiative', target="_blank")
            ],
            nav=True,
            in_navbar=True,
            label="Information",
        ),
    ],
        
    )
)'''

app.layout = dbc.Container(
    [navbar, dash.page_container, dcc.Location(id="url", refresh=True)],
    fluid=True,
)



#Run the dashboard
if __name__ == "__main__":
    app.run_server(debug=True, port=1099)