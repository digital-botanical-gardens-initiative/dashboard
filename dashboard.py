import dash
from dash import dcc
import dash_bootstrap_components as dbc
import webbrowser


app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.config.suppress_callback_exceptions = True

navbar = dbc.NavbarSimple(
    [
        dbc.Button("Home", href="/", color="secondary", className="me-1"),
        dbc.Button("Explore", href="/explore", color="secondary"),
    ],
    brand="Digital Botanical Garden Initiative",
    color="primary",
    dark=True,
    className="mb-2",
)

app.layout = dbc.Container(
    [navbar, dash.page_container, dcc.Location(id="url", refresh=True)],
    fluid=True,
)



#Run the dashboard
if __name__ == "__main__":
    webbrowser.open_new('http://127.0.0.1:8050')
    app.run_server(debug=True, port=8050)