// exploreRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');
const fetch = require('node-fetch');
const { spawn } = require('child_process');


router.use(express.urlencoded({ extended: true }));


const taxonomyColumns = [
  'organism_taxonomy_01domain',
  'organism_taxonomy_02kingdom',
  'organism_taxonomy_03phylum',
  'organism_taxonomy_04class',
  'organism_taxonomy_05order',
  'organism_taxonomy_06family',
  'organism_taxonomy_07tribe',
  'organism_taxonomy_08genus',
  'organism_taxonomy_09species',
  'organism_taxonomy_10varietas'
];


const getTableColumns = async () => {
  const tableColumns = await db.query(`
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'data'
  `);
  return tableColumns.rows;
};

async function sendDataPython(script, query, smiles, tanimoto = 1) {


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

      resolve(matchingSmiles);
    });
    
    process.on('error', (error) => {
      console.error('Error in Python script:', error);
    });
    
  });
}




async function handleExactMatch(radio, display, smiles, group, max, tanimoto) {

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

  InChi = '{' + await getInChi + '}';
  smiles_arr = '{' + smiles + '}';

  if (radio === 'structure_inchi' || radio === 'structure_smiles') {
    let structure = radio === 'structure_inchi' ? InChi : smiles_arr;

    if (display === 'table') {
      return await getTableData(radio, structure, group, max);
    } else if (display === 'graph') {
      return await getGraphData(radio, structure, group, max);
    }
  }

  return { result: [] };
}



// Function to handle 'sub_search' tab
async function handleSubSearch(radio, display, smiles, group, max, tanimoto) {
  const queryGroup = group ? `WHERE '${group}' = ANY (array[${taxonomyColumns.join(', ')}])` : '';
  const query = `SELECT DISTINCT structure_smiles FROM data ${queryGroup}`;

  const matchingSmiles = await sendDataPython('public/python/substructureSearch.py', query, smiles).catch((error) => {
    console.error("Error while calling Python script:", error);
    });
  const column = 'structure_smiles';

    switch(display) {
      case 'table':
        return await getTableData(column, matchingSmiles, group, max);
      case 'graph':
        return await getGraphData(column, matchingSmiles, group, max);
      default:
        throw new Error(`Unknown display type: ${display}`);
    }
  }


// Function to handle 'sim_search' tab

async function handleSimSearch(radio, display, smiles, group, max, tanimoto) {
  const queryGroup = group ? `WHERE '${group}' = ANY (array[${taxonomyColumns.join(', ')}])` : '';
  const query = `SELECT DISTINCT structure_smiles FROM data ${queryGroup}`;

  const matchingSmiles = await sendDataPython('public/python/tanimotoSearch.py',query, smiles, tanimoto).catch((error) => {
    console.error("Error while calling Python script:", error);
});
  const column = 'structure_smiles';

switch(display) {
  case 'table':
    return await getTableData(column, matchingSmiles, group, max);
  case 'graph':
    return await getGraphData(column, matchingSmiles, group, max);
  default:
    throw new Error(`Unknown display type: ${display}`);
}
}

async function getTableData(column, structure, group, max) {
  let query = `SELECT * FROM data WHERE ${column} = ANY($1)`;
  let params = [structure];

  if (group) {
    query += ` AND $2 = ANY (array[${taxonomyColumns.join(', ')}])`;
    params.push(group);
  }

  query += ' LIMIT $' + (params.length + 1);
  params.push(max);

  const result = await db.query(query, params);
  return { result };
}

async function getGraphData(column, structure, group, max) {
  let query = `SELECT organism_taxonomy_01domain as domain,
           organism_taxonomy_02kingdom as kingdom,
           organism_taxonomy_03phylum as phylum,
           organism_taxonomy_04class as class,
           organism_taxonomy_06family as family,
           organism_taxonomy_07tribe as tribe,
           organism_taxonomy_08genus as genus
           FROM data WHERE ${column} = ANY($1)`;

  let params = [structure];

  if (group) {
    query += ` AND $2 = ANY (array[${taxonomyColumns.join(', ')}])`;
    params.push(group);
  }

  query += ' LIMIT $' + (params.length + 1);
  params.push(max);

  const graphData = await db.query(query, params);
  const totalCount = graphData.rowCount;
  const result = buildGraphData(graphData.rows);

  return { result, totalCount };
}


function buildGraphData(rows) {
  const root = { name: "Root", children: [] };

  rows.forEach(row => {
    let currentLevel = root.children;
    [row.domain, row.kingdom, row.phylum, row.class, row.family, row.tribe, row.genus].forEach((level, i, arr) => {
      let existingPath = currentLevel.find(d => d.name === level);
      if (existingPath) {
        currentLevel = existingPath.children;
      } else {
        const newPath = { name: level, children: [] };
        if (i === arr.length - 1) {
          // This is a leaf node, so give it a size.
          newPath.size = 1;
        }
        currentLevel.push(newPath);
        currentLevel = newPath.children;
      }
    });
  });

  return root;
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
    let display = 'table'

    if (req.method === 'POST') {
      const smiles = req.body.smiles;
      display = req.body.display;
      const max = req.body.maxNum;
      const group = req.body.taxo;
      const activeTab = req.body.activeTab;
      const radio = req.body.radio;
      const tanimoto = req.body.tanimoto / 100 ;

      // Get the handler for the active tab
      const tabHandler = tabHandlers[activeTab];
      
      if (tabHandler) {
      
        const results = await tabHandler(radio, display, smiles, group, max, tanimoto) ;
        console.log(results);

        if (display === 'table'){
          console.log(columns);
          res.render('exploreStructure', { columns, results: results.result.rows, hits: results.result.rows.length , display: display});
        } else if (display === 'graph'){
          hits = results.totalCount;
          console.log(hits);
          res.render('exploreStructure', { columns, results: results.result, hits: hits  , display: display});
        }
      } else {
        // Handle unknown tab
        res.send('Unknown tab');
      }
    } else {
      res.render('exploreStructure', { columns, hits: 0 , display: display});
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