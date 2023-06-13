// exploreRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');
const fetch = require('node-fetch');
const { spawn } = require('child_process');
/* const redis = require('redis');
const { promisify } = require('util');

const client = redis.createClient(); // This line may need to be configured based on your Redis setup

// Convert Redis client methods to Promises, to allow us to use async/await syntax
client.getAsync = promisify(client.get).bind(client);
client.setAsync = promisify(client.set).bind(client); */

router.use(express.urlencoded({ extended: true }));

/* const { newEngine } = require('@comunica/actor-init-sparql');
const { SparqlEndpointFetcher } = require('fetch-sparql-endpoint');
const fetcher = new SparqlEndpointFetcher();
const myEngine = newEngine(); */


const getTableColumns = async () => {
  const tableColumns = await db.query(`
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'data'
  `);
  return tableColumns.rows;
};

async function sendDataPython(script, query, smiles, tanimoto = 1) {
/*   console.time('sendDataPython');
  let cacheKey = `sendDataPython:${query}:${smiles}`;
  let cachedResult = await client.getAsync(cacheKey);
  
  if (cachedResult) {
    console.log("Cache hit for", cacheKey);
    return JSON.parse(cachedResult);
  }

  console.log("Cache miss for", cacheKey); */

  return new Promise(async (resolve, reject) => {
    console.time('db.query');
    // Fetch SMILES strings from the database
    const smilesList = await db.query(query);
    console.timeEnd('db.query');
    
    const smilesListMapped = smilesList.rows.map(row => row.structure_smiles);
    
    console.time('spawn process');
    // Call Python script to perform substructure search
    const process = spawn('python3', [script, smiles, tanimoto]);
    process.stdin.write(smilesListMapped.join('\n'));
    process.stdin.on('finish', () => {
      console.log('Finished writing to stdin');
    });
    process.stdin.end();
    
    let matchingSmiles = '';
    process.stdout.on('data', (data) => {
      matchingSmiles += data.toString();
    });
    
    process.stderr.on('data', (data) => {
      console.error(`Python stderr: ${data}`);
    });
    
    process.stdout.on('end', () => {
      // Split the output into a list of SMILES strings
      matchingSmiles = matchingSmiles.split('\n');
      console.log('Python process finished');
      
      // Cache the result before resolving
      /* client.setAsync(cacheKey, JSON.stringify(matchingSmiles))
        .catch(err => console.error('Error caching data:', err));
       */
      resolve(matchingSmiles);
    });
    
    process.on('error', (error) => {
      console.error('Error in Python script:', error);
    });
    
  });
}




async function handleExactMatch(radio, smiles, group, max, tanimoto) {

  let InChi;
  // Promisify the Python spawn process
  const getInChi = new Promise((resolve, reject) => {
    // Call Python script to convert SMILES to InChI
    const process = spawn('python3', ['public/python/smileToInchi.py', smiles]);

    process.stdout.on('data', (data) => {
      // 'data' is the InChI string returned from Python script
      InChi = data.toString();
    });

    process.on('close', (code) => {
      if (code !== 0) {
        reject(new Error(`Python script exited with code ${code}`));
      } else {
        resolve(InChi);
      }
    });
  });

  InChi = await getInChi;
  
  if (radio === 'structure_inchi' || radio === 'structure_smiles') {
    let structure = radio === 'structure_inchi' ? InChi : smiles;
    if (group) {
      let result = await db.query(`
      SELECT * FROM data
      WHERE ${radio} = $1 AND 
      $2 = ANY (array[organism_taxonomy_01domain, 
                      organism_taxonomy_02kingdom, 
                      organism_taxonomy_03phylum,
                      organism_taxonomy_04class,
                      organism_taxonomy_05order,
                      organism_taxonomy_06family,
                      organism_taxonomy_07tribe,
                      organism_taxonomy_08genus,
                      organism_taxonomy_09species,
                      organism_taxonomy_10varietas])
      LIMIT $3
      `, [structure, group, max]);

      return result;
    } else {
      let result = await db.query(`SELECT * FROM data WHERE ${radio} = $1 LIMIT $2`, [structure, max]);
  
      /* if (result) {
        // If the database query is successful, cache the result
        await client.setAsync(cacheKey, JSON.stringify(result));
      } */
      
      return result;    
    } 
  }
  return { rows: [] };
}


// Function to handle 'sub_search' tab
async function handleSubSearch(radio, smiles, group, max, tanimoto) {
/*   let cacheKey = `subSearch:${radio}:${smiles}:${group}:${max}`;
  let cachedResult = await client.getAsync(cacheKey);

  if (cachedResult) {
    console.log("Cache hit for", cacheKey);
    return JSON.parse(cachedResult);
  }
  
  console.log("Cache miss for", cacheKey); */

  if (group) {
    const query = `SELECT structure_smiles FROM data
                                    WHERE ${group} = ANY (array[organism_taxonomy_01domain, 
                                                    organism_taxonomy_02kingdom, 
                                                    organism_taxonomy_03phylum,
                                                    organism_taxonomy_04class,
                                                    organism_taxonomy_05order,
                                                    organism_taxonomy_06family,
                                                    organism_taxonomy_07tribe,
                                                    organism_taxonomy_08genus,
                                                    organism_taxonomy_09species,
                                                    organism_taxonomy_10varietas])`;

    const matchingSmiles = await sendDataPython('public/python/substructureSearch.py',query, smiles).catch((error) => {
      console.error("Error while calling Python script:", error);
      });

    return await db.query(`SELECT structure_smiles FROM data
                                    WHERE structure_smiles = ANY($1) AND
                                    $2 = ANY (array[organism_taxonomy_01domain, 
                                                    organism_taxonomy_02kingdom, 
                                                    organism_taxonomy_03phylum,
                                                    organism_taxonomy_04class,
                                                    organism_taxonomy_05order,
                                                    organism_taxonomy_06family,
                                                    organism_taxonomy_07tribe,
                                                    organism_taxonomy_08genus,
                                                    organism_taxonomy_09species,
                                                    organism_taxonomy_10varietas])`, [matchingSmiles, group, max]);
  } else {
    const query = `SELECT structure_smiles FROM data`;
    const matchingSmiles = await sendDataPython('public/python/substructureSearch.py',query, smiles).catch((error) => {
      console.error("Error while calling Python script:", error);
    });

    return await db.query(`SELECT * FROM data WHERE structure_smiles = ANY($1) LIMIT $2`, [matchingSmiles, max]);
  }
}


// Function to handle 'sim_search' tab
async function handleSimSearch(radio, smiles, group, max, tanimoto) {

  if (group) {
    const query = `SELECT DISTINCT structure_smiles FROM data
                                    WHERE '${group}' = ANY (array[organism_taxonomy_01domain, 
                                                    organism_taxonomy_02kingdom, 
                                                    organism_taxonomy_03phylum,
                                                    organism_taxonomy_04class,
                                                    organism_taxonomy_05order,
                                                    organism_taxonomy_06family,
                                                    organism_taxonomy_07tribe,
                                                    organism_taxonomy_08genus,
                                                    organism_taxonomy_09species,
                                                    organism_taxonomy_10varietas])`;

    const matchingSmiles = await sendDataPython('public/python/tanimotoSearch.py',query, smiles, tanimoto).catch((error) => {
      console.error("Error while calling Python script:", error);
      });

    return await db.query(`SELECT * FROM data
                                    WHERE structure_smiles = ANY($1) AND
                                    $2 = ANY (array[organism_taxonomy_01domain, 
                                                    organism_taxonomy_02kingdom, 
                                                    organism_taxonomy_03phylum,
                                                    organism_taxonomy_04class,
                                                    organism_taxonomy_05order,
                                                    organism_taxonomy_06family,
                                                    organism_taxonomy_07tribe,
                                                    organism_taxonomy_08genus,
                                                    organism_taxonomy_09species,
                                                    organism_taxonomy_10varietas])
                                    LIMIT $3`
                                    , [matchingSmiles, group, max]);
  } else {
    const query = `SELECT DISTINCT structure_smiles FROM data`;
    const matchingSmiles = await sendDataPython('public/python/tanimotoSearch.py',query, smiles, tanimoto).catch((error) => {
      console.error("Error while calling Python script:", error);
    });

    return await db.query(`SELECT * FROM data WHERE structure_smiles = ANY($1) LIMIT $2`, [matchingSmiles, max]);
  }
}
// Add more handler functions for other tabs...

// A mapping from tabs to their handler functions
const tabHandlers = {
  'exact_match': handleExactMatch,
  'sub_search': handleSubSearch,
  'sim_search': handleSimSearch
  // Add more mappings for other tabs...
};


router.get('/explore', async (req, res) => {
  try {
    const columns = await getTableColumns();
    res.render('explore', { columns });
  } catch (err) {
    console.error(err);
    res.send('Error while fetching columns names');
  }
});

router.all('/explore/text', async (req, res) => {
  try {
    const columns = await getTableColumns();

    if (req.method === 'POST') {
      const column = req.body.column;
      const searchTerm = req.body.search;
      let results = {}
      if (req.body.exact === 'true') {
        results = await db.query(`SELECT * FROM data WHERE ${column} = $1`, [searchTerm]);
      } else {
        results = await db.query(`SELECT * FROM data WHERE ${column} ILIKE $1`, [`%${searchTerm}%`]);
      }
      
      res.render('exploreText', { columns, results: results.rows , hits: results.rows.length}); 
    } else {
      res.render('exploreText', { columns , hits: 0});
    }
  } catch (err) {
    console.error(err);
    res.send('Error while fetching columns names');
  }
});


router.all('/explore/structure', async (req, res) => {
  try {
    const columns = await getTableColumns();

    if (req.method === 'POST') {
      const smiles = req.body.smiles;
      const max = req.body.maxNum;
      const group = req.body.taxo;
      const activeTab = req.body.activeTab;
      const radio = req.body.radio;
      const tanimoto = req.body.tanimoto / 100 ;

      // Get the handler for the active tab
      const tabHandler = tabHandlers[activeTab];
      
      if (tabHandler) {
      
        const results = await tabHandler(radio, smiles, group, max, tanimoto);

        res.render('exploreStructure', { columns, results: results.rows, hits: results.rows.length });

        
      } else {
        // Handle unknown tab
        res.send('Unknown tab');
      }
    } else {
      res.render('exploreStructure', { columns, hits: 0 });
    }
  } catch (err) {
    console.error(err);
    res.send('Oops... something went wrong!');
  }
});


router.all('/explore/SPARQL', async (req,res) => {
  try {
    res.render('exploreSparql');
  } catch (err) {
    console.error(err);
    res.send('Oops... something went wrong!');
  }
})

/* router.post('/sparql', async (req, res) => {
  const query = req.body.query;
  
  try {
    const result = await myEngine.query(query, {
      sources: [{ type: 'sparql', value: endpointUrl }],
    });

    const bindings = await result.bindings();
    const results = bindings.map(binding => binding.toObject());
    res.render('sparql', { results });
  } catch (err) {
    console.error(err);
    res.send('Error while executing SPARQL query');
  }
}); */

module.exports = router;