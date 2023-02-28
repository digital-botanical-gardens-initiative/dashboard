import dash
import dash_labs as dl
import dash_bootstrap_components as dbc


icon = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.3/font/bootstrap-icons.css"

# meta_tags are required for the app layout to be mobile responsive
app = dash.Dash(
    __name__,use_pages=True, pages_folder='', external_stylesheets=[dbc.themes.BOOTSTRAP,icon]
)

server = app.server