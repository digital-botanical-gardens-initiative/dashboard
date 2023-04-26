import os
from sqlalchemy import create_engine
from rdkit import Chem
from rdkit.Chem import AllChem
from sqlalchemy.orm import sessionmaker
from sqlalchemy import insert
from dask.delayed import delayed
from sqlalchemy import text
import pandas as pd
import dask.dataframe as dd
import numpy as np
import dask.array as da
import psutil
from sklearn.manifold import TSNE
import sqlalchemy


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

def print_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    rss = mem_info.rss / (1024 ** 2)  # Convert to MB
    print(f"Memory usage: {rss:.2f} MB")


# Define the table
umap_table = sqlalchemy.Table(
    'umap',
    sqlalchemy.MetaData(),
    sqlalchemy.Column('index', sqlalchemy.Integer, primary_key=True),
    sqlalchemy.Column('UMAP1', sqlalchemy.Float),
    sqlalchemy.Column('UMAP2', sqlalchemy.Float),
    sqlalchemy.Column('class_name', sqlalchemy.String)
)

db_connection_str = f'postgresql://{db_user}:{db_pwd}@localhost/{db_name}'
engine = create_engine(db_connection_str)
Session = sessionmaker(bind=engine)

query_in = f"SELECT structure_smiles, structure_taxonomy_npclassifier_03class FROM {db_name}"

print("Before loading data:")
print_memory_usage()

# Define a function to read chunks of data from the database
@delayed
def load_chunk(chunk_id, chunk_size):
    offset = chunk_id * chunk_size
    chunk_query = f"{query_in} LIMIT {chunk_size} OFFSET {offset}"
    with engine.connect() as conn:
        chunk = pd.read_sql(text(chunk_query), con=conn)
        chunk['ECFP'] = chunk['structure_smiles'].apply(lambda x: np.array(smiles_to_ecfp(x), dtype=np.float64))
        return chunk

# Set the number of chunks and chunk size based on your dataset size and available memory
num_chunks = 5
chunk_size = 20000

# Read data from the PostgreSQL table into Dask DataFrame using the delayed function
chunks = [load_chunk(chunk_id, chunk_size) for chunk_id in range(num_chunks)]
ddf = dd.from_delayed(chunks)

print("After processing data:")
print_memory_usage()

#ddf['ECFP'] = ddf['structure_smiles'].apply(lambda x: np.array(smiles_to_ecfp(x), dtype=np.float64), meta=(None, np.float64))
ddf = ddf.persist()  # Persist the ECFP column in memory

print('apply done')

# Convert the ECFP column to a numerical feature matrix
ecfp_series = ddf['ECFP'].compute()
X = np.stack(ecfp_series.values)

print("Before applying t-SNE:")
print_memory_usage()

# Apply t-SNE dimensionality reduction
tsne = TSNE(n_components=2)
tsne_embedding = tsne.fit_transform(X)

print("After applying t-SNE:")
print_memory_usage()

# Create a DataFrame for the t-SNE embedding with the family labels
class_labels = ddf['structure_taxonomy_npclassifier_03class'].compute().reset_index(drop=True)
tsne_df = pd.DataFrame(tsne_embedding, columns=['UMAP1', 'UMAP2'])
tsne_df['class_name'] = class_labels

conn = engine.connect()
trans = conn.begin()  # Start a new transaction

print('connection opened')

# Write each row of the DataFrame to the database
for idx, row in tsne_df.iterrows():
    stmt = insert(umap_table).values(
        index=idx,
        UMAP1=row['UMAP1'],
        UMAP2=row['UMAP2'],
        class_name=row['class_name']
    )
    conn.execute(stmt)

trans.commit()  # Commit the transaction
conn.close()

print('connection closed')


with engine.connect() as conn:
    result = conn.execute(sqlalchemy.text("SELECT COUNT(*) FROM umap"))
    count = result.scalar()
    print(f"Total number of rows in the 'umap' table: {count}")
