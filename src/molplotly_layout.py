from dash import Dash, html, dcc, Input, Output
import plotly.express as px
import pandas as pd
import vaex


# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
import molplotly

# load a DataFrame with smiles
df = pd.read_csv("/home/mwannier/dbgi_dashboard/data/230106_frozen_metadata.csv",nrows=1000)
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

layout = jupyterFig.layout

