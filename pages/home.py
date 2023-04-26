import dash
from dash import html, register_page, dcc
import dash_bootstrap_components as dbc
import plotly.express as px
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from sqlalchemy import create_engine
from dask.delayed import delayed
import pandas as pd
from sqlalchemy.orm import sessionmaker
from sqlalchemy import text
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


# import env variable
from dotenv import load_dotenv
load_dotenv()

db_name=os.getenv('DB_NAME')
db_pwd=os.getenv('DB_PWD')
db_user=os.getenv('DB_USR')


def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    ecfp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return ecfp


# Replace the following with your own database connection information
db_connection_str = f'postgresql://{db_user}:{db_pwd}@localhost/{db_name}'
engine = create_engine(db_connection_str)
Session = sessionmaker(bind=engine)

# Replace 'your_table_name' with the actual table name in your database
query_lotus = f"SELECT * FROM {db_name}"
query_umap = 'SELECT DISTINCT * FROM umap'

# Define a function to read chunks of data from the database
@delayed
def load_chunk(chunk_id, chunk_size, query):
    offset = chunk_id * chunk_size
    chunk_query = f"{query} LIMIT {chunk_size} OFFSET {offset}"
    with engine.connect() as conn:
        return pd.read_sql(text(chunk_query), con=conn)
    

# Set the number of chunks and chunk size based on your dataset size and available memory
num_chunks = 5
chunk_size = 20000

# Read data from the PostgreSQL table into Dask DataFrame using the delayed function
chunks_lotus = [load_chunk(chunk_id, chunk_size, query_lotus) for chunk_id in range(num_chunks)]
ddf = dd.from_delayed(chunks_lotus)

# Read data from the PostgreSQL table into Dask DataFrame using the delayed function
chunks_umap = [load_chunk(chunk_id, chunk_size, query_umap) for chunk_id in range(num_chunks)]
umap_df = dd.from_delayed(chunks_umap)

register_page(__name__, path="/")

fields=['structure_smiles',
        'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 
        'organism_taxonomy_03phylum', 'organism_taxonomy_04class',
        'organism_taxonomy_05order', 'organism_taxonomy_06family',
        'organism_taxonomy_07tribe','organism_taxonomy_08genus',
        'structure_taxonomy_npclassifier_01pathway','structure_taxonomy_npclassifier_02superclass',
        'structure_taxonomy_npclassifier_03class','structure_nameTraditional']

df_orga = ddf.query('organism_taxonomy_01domain == organism_taxonomy_01domain').fillna('NA')
#df_chemo = ddf.query('structure_taxonomy_npclassifier_01pathway == structure_taxonomy_npclassifier_01pathway').fillna('NA')

# Aggregate data for both DataFrames by count
df_orga_agg = df_orga.groupby(fields[1:8]).size().reset_index().rename(columns={0: 'count'})
#df_chemo_agg = df_chemo.groupby(fields[8:]).size().reset_index().rename(columns={0: 'count'})


with ProgressBar():
    #fig_chemo = px.treemap(df_chemo_agg.compute(), path=fields[8:], values='count')
    fig_orga = px.treemap(df_orga_agg.compute(), path=fields[1:8], values='count')
    fig_chemo = px.scatter(umap_df, x=0, y=1, color=2, title='UMAP projection of the molecular dataset', template='simple_white')


layout = html.Div([
                    html.Div([
                            html.H1('Welcome!'),
                            html.H3('Molecular diversity in the DBGI'),
                            dcc.Graph(figure=fig_chemo),
                            dbc.Button('Browse Molecular diversity',id='browse-mol', href='/element'),
                        ]),
                    html.Br(),
                    html.Div([
                            html.H3("Organism diversity in the DBGI"),
                            dcc.Graph(figure=fig_orga),
                            dbc.Button('Browse Organism diversity',id='browse-orga', href='/organism'),
                        ]),
                    html.Br()
                ])
