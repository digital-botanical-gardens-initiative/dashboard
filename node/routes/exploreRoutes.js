// exploreRoutes.js
// Importing necessary modules
const express = require('express');
const router = express.Router();
const db = require('./db');
const fetch = require('node-fetch');
const { spawn } = require('child_process');
const dotenv = require("dotenv");
const { rawListeners } = require('process');
dotenv.config();

// Importing necessary modules from 'graphdb' package
const {ServerClient, GraphDBServerClient, ServerClientConfig} = require('graphdb').server;
const {RepositoryClientConfig, RDFRepositoryClient} = require('graphdb').repository;
const {RDFMimeType} = require('graphdb').http;
const {SparqlJsonResultParser} = require('graphdb').parser;
const {GetQueryPayload, QueryType} = require('graphdb').query;

// Middleware to parse URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

// Defining the columns related to taxonomy
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

// Function to retrieve column names from the 'data' table in the 'public' schema
const getTableColumns = async () => {
  const tableColumns = await db.query(`
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'data'
  `);
  return tableColumns.rows;
};

// Function to execute a Python script with some inputs and return the results
async function sendDataPython(script, query, smiles, tanimoto = 1) {
  return new Promise(async (resolve, reject) => {
    // Timing the db query execution
    console.time('db.query');
    const smilesList = await db.query(query);  // Fetch SMILES strings from the database
    console.timeEnd('db.query');
    
    // Mapping through the response to extract only the SMILES structure
    const smilesListMapped = smilesList.rows.map(row => row.structure_smiles);
    
    // Calling a Python script to perform a substructure search
    console.time('spawn process');
    const process = spawn('python3', [script, smiles, tanimoto]);
    process.stdin.write(smilesListMapped.join('\n'));
    process.stdin.on('finish', () => {
      console.log('Finished writing to stdin');
    });
    process.stdin.end();
    
    let matchingSmiles = '';
    process.stdout.on('data', (data) => {
      // Appending the output data to 'matchingSmiles'
      matchingSmiles += data.toString();
    });
    
    // Logging any errors in the Python script execution
    process.stderr.on('data', (data) => {
      console.error(`Python stderr: ${data}`);
    });
    
    process.stdout.on('end', () => {
      // Splitting the output into a list of SMILES strings
      matchingSmiles = matchingSmiles.split('\n');
      console.log('Python process finished');
      resolve(matchingSmiles);  // Resolve the promise
    });
    
    process.on('error', (error) => {
      console.error('Error in Python script:', error);
    });
  });
}

// Function to handle 'exact_search' tab
async function handleExactMatch(radio, display, smiles, group, max, tanimoto) {
  // Convert SMILES to InChI using Python
  let InChi;
  const getInChi = new Promise((resolve, reject) => {
    const process = spawn('python3', ['public/python/smileToInchi.py', smiles]);
    process.stdout.on('data', (data) => {
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


// Function for finding exact text matches in the database
async function exactTextMatch(display, column, text, max){
  var group = ''

  // Format text as string array for PostgreSQL
  var text = '{' + text + '}';

  // Depending on the display type, return table or graph data
  switch(display) {
    case 'table':
      return await getTableData(column, text, group, max);
    case 'graph':
      return await getGraphData(column, text, group, max);
    default:
      throw new Error(`Unknown display type: ${display}`);
  }
}

async function handleSPARQLExactMatch(radio, display, smiles, group, max, tanimoto){
    // Convert SMILES to InChI using Python
    let InChi;
    const getInChi = new Promise((resolve, reject) => {
      const process = spawn('python3', ['public/python/smileToInchiKey.py', smiles]);
      process.stdout.on('data', (data) => {
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
       let query = `PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
                    select ?InchiKey ?smiles ?NPCpathway ?NPCsuperclass ?NPCclass ?wd_id `;

        if(radio === 'structure_inchi'){
          query += `where {
            <https://enpkg.commons-lab.org/kg/${structure}> enpkg:has_smiles ?smiles ;
                                                              enpkg:has_wd_id ?wd_id ;
                                                              enpkg:has_npc_pathway ?NPCpathway ;
                                                              enpkg:has_npc_superclass ?NPCsuperclass ;
                                                              enpkg:has_npc_class ?NPCclass .
            BIND(<https://enpkg.commons-lab.org/kg/${structure}> AS ?InchiKey) }
            LIMIT ${max}`;
        } else if(radio === 'structure_smiles'){
          query += `where {
                  ?InchiKey enpkg:has_smiles "${structure}" ;
                                                              enpkg:has_wd_id ?wd_id ;
                                                              enpkg:has_npc_pathway ?NPCpathway ;
                                                              enpkg:has_npc_superclass ?NPCsuperclass ;
                                                              enpkg:has_npc_class ?NPCclass .
            BIND("${structure}" AS ?smiles) }
            LIMIT ${max}`;
        }
      // Query the GraphDB using the SPARQL query and get the results
      const graph = await queryGraphDB(query);

      // Render the 'exploreSparql' page with the query results
      return {results: graph.nodes, headers: graph.headers};  
 
    }
    return { result: [] };
}

async function handleSPARQLSubSearch(radio, display, smiles, group, max, tanimoto){
    console.log(smiles);
    let query = `PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
                PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                PREFIX wdt: <http://www.wikidata.org/prop/direct/>
                PREFIX wd: <http://www.wikidata.org/entity/>
                PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/>
                PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
                SELECT DISTINCT ?ik ?smiles ?NPCpathway ?NPCsuperclass ?NPCclass ?wd_id
                WHERE
                { ?ik enpkg:has_smiles ?smiles ;
                  enpkg:has_npc_pathway ?NPCpathway ;
                  enpkg:has_npc_superclass ?NPCsuperclass ;
                  enpkg:has_npc_class ?NPCclass ;
                  enpkg:has_wd_id ?wd_id .
                  SERVICE idsm:wikidata {
                    VALUES ?SUBSTRUCTURE {
                        "${smiles}" # Structure given by the user
                    }
                    ?wd_id sachem:substructureSearch _:b16.
                    _:b16 sachem:query ?SUBSTRUCTURE.
                  }      
                }
                LIMIT ${max}`;

        // Query the GraphDB using the SPARQL query and get the results
        const graph = await queryGraphDB(query);
        console.log(graph);

        // Render the 'exploreSparql' page with the query results
        return {results: graph.nodes, headers: graph.headers};  
}

async function handleSPARQLSimSearch(radio, display, smiles, group, max, tanimoto){
      let query = `PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
      PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
      PREFIX wdt: <http://www.wikidata.org/prop/direct/>
      PREFIX wd: <http://www.wikidata.org/entity/>
      PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/>
      PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
      SELECT DISTINCT ?ik ?smiles ?NPCpathway ?NPCsuperclass ?NPCclass ?wd_id
      WHERE
      { 
          ?extract rdf:type enpkg:LabExtract .
              ?extract enpkg:has_LCMS ?lcms .
                  ?lcms enpkg:has_lcms_feature_list ?feature_list .
                  ?feature_list enpkg:has_lcms_feature ?feature .
                      ?feature enpkg:has_sirius_annotation|enpkg:has_isdb_annotation ?annotation . 
                      ?annotation enpkg:has_InChIkey2D ?ik2d .
                          ?ik2d enpkg:has_smiles ?smiles .
                          ?ik2d enpkg:is_InChIkey2D_of ?ik .
                            ?ik enpkg:has_smiles ?smiles ;
                                enpkg:has_npc_pathway ?NPCpathway ;
                                enpkg:has_npc_superclass ?NPCsuperclass ;
                                enpkg:has_npc_class ?NPCclass ;
                                enpkg:has_wd_id ?wd_id .
                              SERVICE idsm:wikidata {
                                ?smiles sachem:similarCompoundSearch [
                                  sachem:query "${smiles}";
                                  sachem:cutoff "${tanimoto}"^^xsd:double ].
                              }      
      }
      LIMIT ${max}`;
      console.log(query);
      // Query the GraphDB using the SPARQL query and get the results
      const graph = await queryGraphDB(query);

      // Render the 'exploreSparql' page with the query results
      return {results: graph.nodes, headers: graph.headers};  
}

// Function to fetch data for table display
async function getTableData(column, structure, group, max) {
  // SQL query to select data where column value is in the provided structure
  let query = `SELECT * FROM data WHERE ${column} = ANY($1)`;
  let params = [structure];

  // If group is provided, append to SQL query
  if (group) {
    query += ` AND $2 = ANY (array[${taxonomyColumns.join(', ')}])`;
    params.push(group);
  }

  // Limit the query results
  query += ' LIMIT $' + (params.length + 1);
  params.push(max);

  // Execute the query
  const result = await db.query(query, params);
  return { result };
}


// Function to fetch data for graph display
async function getGraphData(column, structure, group, max) {
  // SQL query to select taxonomical data where column value is in the provided structure
  let query = `SELECT organism_taxonomy_01domain as domain,
           organism_taxonomy_02kingdom as kingdom,
           organism_taxonomy_03phylum as phylum,
           organism_taxonomy_04class as class,
           organism_taxonomy_06family as family,
           organism_taxonomy_07tribe as tribe,
           organism_taxonomy_08genus as genus,
           organism_taxonomy_09species as species,
           organism_taxonomy_10varietas as varietas
           FROM data WHERE ${column} = ANY($1)`;

  let params = [structure];

  // If group is provided, append to SQL query
  if (group) {
    query += ` AND $2 = ANY (array[${taxonomyColumns.join(', ')}])`;
    params.push(group);
  }

  // Limit the query results
  query += ' LIMIT $' + (params.length + 1);
  params.push(max);

  // Execute the query
  const graphData = await db.query(query, params);
  const totalCount = graphData.rowCount;

  // Build graph data from the results
  const result = buildGraphData(graphData.rows);

  return { result, totalCount };
}

// Function to convert rows of taxonomy data into hierarchical structure for graph visualization
function buildGraphData(rows) {
  // Initialize root of the graph
  const root = { name: "Root", children: [], size: 0 };

  // For each row of data
  rows.forEach(row => {
    let currentLevel = root.children;
    // For each level of taxonomy, add a node to the graph
    [row.domain, row.kingdom, row.phylum, row.class, row.family, row.tribe, row.genus, row.species, row.varietas].forEach((level, i, arr) => {
      // Look for an existing node at the current level
      let existingPath = currentLevel.find(d => d.name === level);
      if (existingPath) {
        // If found, descend to its children and increment size
        existingPath.size++;
        currentLevel = existingPath.children;
      } else {
        // If not found, create a new node and descend to it
        const newPath = { name: level, children: [], size: 1 };
        currentLevel.push(newPath);
        currentLevel = newPath.children;
      }
    });
  });

  return root;
}


// Function to query a GraphDB and return the data
async function queryGraphDB(query) {
  // Define the endpoint, readTimeout and writeTimeout for GraphDB
  const endpoint = process.env.ENDPOINT_GRAPHDB;
  const readTimeout = 300000;
  const writeTimeout = 300000;
  
  // Set up the configuration for the RDFRepositoryClient
  const config = new RepositoryClientConfig(endpoint)
      .setEndpoints([process.env.REPO_GRAPHDB])
      .setHeaders({
        'Accept': RDFMimeType.TURTLE
      })
      .setReadTimeout(readTimeout)
      .setWriteTimeout(writeTimeout);
      
  // Create a new RDFRepositoryClient with the config
  const repository = new RDFRepositoryClient(config);

  // Register a new SparqlJsonResultParser
  repository.registerParser(new SparqlJsonResultParser());

  // Define the payload for the query
  const payload = new GetQueryPayload()
                        .setQuery(query)
                        .setQueryType(QueryType.SELECT)
                        .setResponseType(RDFMimeType.SPARQL_RESULTS_JSON)
                        .setLimit(100);

  // Return a promise that resolves with the data from the GraphDB
  return new Promise((resolve, reject) => {
    const nodes = [];
    let headers = [];

    // Query the repository with the payload
    repository.query(payload).then((stream) => {
      stream.on('data', (bindings) => {
        const node = {};
        // For each binding in the data, add it to the node
        for (const binding in bindings) {
          if (bindings.hasOwnProperty(binding)) {
            node[binding] = bindings[binding].value;  // use binding as id and bindings[binding].value as label
            if (!headers.includes(binding)) {
              headers.push(binding);  // add the binding to the headers if it's not already there
            }
          }
        }
        nodes.push(node);
      });
      stream.on('end', () => {
        // Resolve the promise with the graph data and the headers when the stream ends
        resolve({nodes: nodes, headers: headers});
      });

      stream.on('error', (error) => {
        // Reject the promise if there's an error
        reject(error);
      });
    });
  });
}

// Handlers mapping object for each tab in the application UI
const tabHandlers = {
  'exact_match': [handleExactMatch, handleSPARQLExactMatch],
  'sub_search': [handleSubSearch, handleSPARQLSubSearch],
  'sim_search': [handleSimSearch, handleSPARQLSimSearch]
  // Add more mappings for other tabs...
};

// Router for '/explore' endpoint, renders the explore page with table columns 
router.get('/explore', async (req, res) => {
  try {
    const columns = await getTableColumns();
    res.render('explore', { columns });
  } catch (err) {
    console.error(err);
    res.send('Error while fetching columns names');
  }
});

// Router for '/explore/text' endpoint, handles both GET and POST requests
router.all('/explore/text', async (req, res) => {
  try {
    const columns = await getTableColumns();
    let display = 'table';

    // If it's a POST request, process the search inputs
    if (req.method === 'POST') {
      const column = req.body.column;
      const searchTerm = req.body.search;
      display = req.body.display;
      const max = req.body.maxNum;
      const datasource = req.body.datasource;

      // Execute the text search and log results
      let results = await exactTextMatch(display, column, searchTerm, max);

      // Render the results in the appropriate format
      if (display === 'table'){
        res.render('exploreText', { columns, results: results.result.rows, hits: results.result.rows.length , display: display});
      } else if (display === 'graph'){
        hits = results.totalCount;
        res.render('exploreText', { columns, results: results.result, hits: hits  , display: display});
      }

    } else {
      // If it's not a POST request, render the page with default settings
      res.render('exploreText', { columns , hits: 0, display});
    }
  } catch (err) {
    // If there's an error, log it and send an error message
    console.error(err);
    res.send('Error while fetching columns names');
  }
});


// Define a router for the '/explore/structure' path which handles all types of HTTP requests
router.all('/explore/structure', async (req, res) => {
  try {
    // Get the table columns
    const columns = await getTableColumns();
    let display = 'table';
    let datasource = 'lotus';

    // Check if the HTTP request is a POST request
    if (req.method === 'POST') {
      // Extract information from the request body
      const smiles = req.body.smiles;
      display = req.body.display;
      const max = req.body.maxNum;
      const group = req.body.taxo;
      const activeTab = req.body.activeTab;
      const radio = req.body.radio;
      const tanimoto = req.body.tanimoto / 100 ;
      datasource = req.body.data_source;

      let tabHandler

      if (datasource === 'lotus'){ 
        // Get the appropriate handler function for the active tab
        tabHandler = tabHandlers[activeTab][0];
      } else if (datasource === 'dbgi'){
        tabHandler = tabHandlers[activeTab][1];
      } else {
        res.send('Unknown source of Data');
      }
      
      if (tabHandler) {
        // If the handler function exists, call it and get the results
        const results = await tabHandler(radio, display, smiles, group, max, tanimoto) ;

        // Based on the 'display' value, render the response
        if (display === 'table'){
          if (datasource === 'lotus'){
            res.render('exploreStructure', { columns, results: results.result.rows, hits: results.result.rows.length , display: display, source: datasource});
          } else if (datasource === 'dbgi') {
            res.render('exploreStructure', {results: results.results, headers: results.headers, hits: 0 , display:display, source: datasource}); 
          }
        } else if (display === 'graph'){
          hits = results.totalCount;
          res.render('exploreStructure', { columns, results: results.result, hits: hits  , display: display});
        }
      } else {
        // If the handler function does not exist, send an error message
        res.send('Unknown tab');
      }
    } else {
      // If the request is not a POST request, render the default 'exploreStructure' page
      res.render('exploreStructure', { columns, hits: 0 , display: display, source: datasource});
    }
  } catch (err) {
    // If there's any error, log it and send an error message
    console.error(err);
    res.send('Oops... something went wrong!');
  }
});

// Define a router for the '/explore/SPARQL' path which handles all types of HTTP requests
router.all('/explore/SPARQL', async (req, res) => {
  try {
    // Check if the HTTP request is a POST request
    if (req.method === 'POST') {
      // Extract the SPARQL query from the request body
      const query = req.body.query;
      
      // Query the GraphDB using the SPARQL query and get the results
      const graph = await queryGraphDB(query);

      // Render the 'exploreSparql' page with the query results
      res.render('exploreSparql', {results: graph.nodes, headers: graph.headers}); 
    } else {
      // If the request is not a POST request, render the default 'exploreSparql' page
      res.render('exploreSparql');
    }
  } catch (error) {
    // If there's any error, log it and send a server error message
    console.error(error);
    res.status(500).send("Server Error");
  }
});

// Export the router module
module.exports = router;
