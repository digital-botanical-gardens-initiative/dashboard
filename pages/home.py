import dash
from dash import html, register_page
import dash_bootstrap_components as dbc


register_page(__name__, path="/")


layout = html.H1('Welcome!')




