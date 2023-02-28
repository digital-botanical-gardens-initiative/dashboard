import dash
from dash import html

header = html.H3('Welcome to home page!')

layout = html.Div([
        header,
    ])

dash.register_page('Home', path='/', layout=layout, icon="bi bi-house")


