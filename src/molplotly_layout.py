import dash
from dash import Dash, html, dcc, Input, Output, no_update
import plotly.express as px
import pandas as pd
from io import BytesIO
import base64

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D

from application import app


# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
import molplotly

# load a DataFrame with smiles
df = pd.read_csv("/home/mwannier/dashboard/data/230106_frozen_metadata.csv")
df['y_pred'] = df['structure_xlogp']
df['y_true'] = df['structure_exact_mass']
#print(df)

# generate a scatter plot
fig = px.scatter(df, x='y_true', y='y_pred')

# add molecules to the plotly graph - returns a Dash app
jupyterFig = molplotly.add_molecules(fig=fig,
                            df=df,
                            smiles_col='structure_smiles_2D',
                            title_col='structure_nameTraditional',
                            )

smiles_col = 'structure_smiles_2D'

if isinstance(smiles_col, str):
    smiles_col = [smiles_col]


layout = html.Div(
        [
            dcc.Graph(id="graph-basic-2", figure=fig, clear_on_unhover=True),   
            dcc.Tooltip(
                id="graph-tooltip", background_color=f"rgba(255,255,255,0.75)"
            ),
        ]
    )

@app.callback(
        output=[
           Output("graph-tooltip", "show"),
           Output("graph-tooltip", "bbox"),
           Output("graph-tooltip", "children"),
        ],
        inputs=[Input("graph-basic-2", "hoverData")],
    )
def display_hover(hoverData):
        mol_alpha = 0.7
        svg_size = 200
        value =  smiles_col
        fontfamily = 'Arial'
        fontsize = 12
        width = 150


        if len(fig.data) != 1:
            colors = {index: x.marker["color"] for index, x in enumerate(fig.data)}
        else:
            colors = {0: "black"}


        if hoverData is None:
            return False, no_update, no_update
        if value is None:
            if smiles_col is not None:
                value = smiles_col
        if isinstance(value, str):
            chosen_smiles = [value]
        else:
            chosen_smiles = value
        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        num = pt["pointNumber"]
        curve_num = pt["curveNumber"]
        if len(fig.data) != 1:
            df_curve = curve_dict[curve_num].reset_index(drop=True)
            df_row = df_curve.iloc[num]
        else:
            df_row = df.iloc[num]
        hoverbox_elements = []

        for col in chosen_smiles:
            smiles = df_row[col]
            buffered = BytesIO()
            if isinstance(smiles, str):
                # Generate 2D SVG if smiles column is a string
                d2d = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
                opts = d2d.drawOptions()
                opts.clearBackground = False
                d2d.DrawMolecule(Chem.MolFromSmiles(smiles))
                d2d.FinishDrawing()
                img_str = d2d.GetDrawingText()
                buffered.write(str.encode(img_str))
                img_str = base64.b64encode(buffered.getvalue())
                img_str = f"data:image/svg+xml;base64,{repr(img_str)[2:-1]}"
            elif isinstance(smiles, Mol):
                # if smiles column is a Mol object, use the 3D coordinates of the mol object
                img = Chem.Draw.MolToImage(smiles)
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue())
                img_str = "data:image/png;base64,{}".format(repr(img_str)[2:-1])
            else:
                raise TypeError(
                    "smiles_col or mol_col not specified with the correct type."
                )
            if len(smiles_col) > 1:
                hoverbox_elements.append(
                    html.H2(
                        f"{col}",
                        style={
                            "color": colors[curve_num],
                            "font-family": fontfamily,
                            "fontSize": fontsize + 2,
                        },
                    )
                )
            hoverbox_elements.append(
                html.Img(
                    src=img_str,
                    style={
                       "width": "100%",
                       "background-color": f"rgba(255,255,255,{mol_alpha})",
                    },
               )
            )

        title_col = 'structure_nameTraditional'
        title = df_row[title_col]
        hoverbox_elements.append(
               html.H4(
                   f"{title}",
                   style={
                       "color": colors[curve_num],
                       "font-family": fontfamily,
                       "fontSize": fontsize,
                   },
               )
           )
        
        children = [
           html.Div(
               hoverbox_elements,
               style={
                   "width": f"{width}px",
                   "white-space": "normal",
               },
           )
       ]

        return True, bbox, children


dash.register_page(' molplotly', path='/molplotly', layout=layout, icon="bi bi-arrow-through-heart-fill")

