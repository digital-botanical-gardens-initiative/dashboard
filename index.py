from dash import html, dcc, Input, Output

# Connect to main app.py file
from app import app, server

# Connect to your app pages
from src import piechart_layout, molplotly_layout


app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div([
        dcc.Link('Page 1|', href='/src/molplotly_layout'),
        dcc.Link('Page 2', href='/src/piechart_layout'),
    ], className="row"),
    html.Div(id='page-content', children=[
        html.H1("Digital Botanical Garden Initiative")
    ])
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/src/molplotly_layout':
        return molplotly_layout.layout
    if pathname == '/src/piechart_layout':
        return piechart_layout.layout



if __name__ == '__main__':
    app.run_server(debug=True)