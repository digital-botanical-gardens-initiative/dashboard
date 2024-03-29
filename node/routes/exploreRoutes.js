/**
 * File/Folder Overview:
 * exploreRoutes.js - Defines routes and utilities for exploring data.
 */

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

// Variable Comments:
// Array defining columns related to taxonomy in the data schema
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

/**
 * Function/Method Comments:
 * An asynchronous function that retrieves column names from the 'data' table in the 'public' schema of the database.
 * 
 * Expected outputs:
 * An array of column names from the specified table and schema.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries. Rethrows errors to allow external error handling.
 */
const getTableColumns = async () => {
  try {
    const tableColumns = await db.query(`
      SELECT column_name
      FROM information_schema.columns
      WHERE table_schema = 'public' AND table_name = 'data'
    `);
    return tableColumns.rows;
  } catch (error) {
    console.error('An error occurred while retrieving table columns:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
};

/**
 * Function/Method Comments:
 * An asynchronous function that executes a Python script to process SMILES strings. The function fetches SMILES data from the database and sends it to the Python script for processing.
 * 
 * Parameters:
 * - script: Path to the Python script to execute.
 * - query: SQL query to fetch SMILES strings from the database.
 * - smiles: The SMILES string to be processed.
 * - tanimoto: A threshold for Tanimoto similarity (default is 1).
 * 
 * Expected outputs:
 * An array of SMILES strings that match or are processed according to the criteria defined in the Python script.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries, Python script execution, and data processing. Rethrows errors to allow external error handling.
 */
async function sendDataPython(script, query, smiles, tanimoto = 1) {
  try {
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
        reject(error); // Reject the promise in case of an error in the python script
      });
    });
  } catch (error) {
    console.error('An error occurred:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function to process and handle exact matches based on the provided parameters. 
 * Primarily converts SMILES to InChI using a Python script.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Data that matches the provided criteria.
 * 
 * Error Handling:
 * Catches and logs errors related to processing and fetching the data. Rethrows errors to allow external error handling.
 */
async function handleExactMatch(radio, display, smiles, group, max, tanimoto) {
  try {
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
  } catch (error) {
    console.error('An error occurred while handling the exact match:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
}

/**
 * Function/Method Comments:
 * An asynchronous function to process and handle substructure searches based on the provided parameters. 
 * Primarily involves querying for data based on the given parameters and matching with SMILES data.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Data that matches the provided criteria.
 * 
 * Error Handling:
 * Catches and logs errors related to processing and fetching the data. Rethrows errors to allow external error handling.
 */
async function handleSubSearch(radio, display, smiles, group, max, tanimoto) {
  try {
    const queryGroup = group ? `WHERE '${group}' = ANY (array[${taxonomyColumns.join(', ')}])` : '';
    const query = `SELECT DISTINCT structure_smiles FROM data ${queryGroup}`;

    const matchingSmiles = await sendDataPython('public/python/substructureSearch.py', query, smiles).catch((error) => {
      console.error("Error while calling Python script:", error);
      // Here you may choose to either rethrow the error or handle it in some other way.
      // For instance, we can return an empty object to signify that no results were found due to an error.
      return {};
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
  } catch (error) {
    console.error('An error occurred while handling the sub search:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function to process and handle similarity searches based on the provided parameters. 
 * Involves querying for data based on the given parameters and comparing SMILES data for similarity.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Data that matches the provided criteria.
 * 
 * Error Handling:
 * Catches and logs errors related to processing and fetching the data. Rethrows errors to allow external error handling.
 */
async function handleSimSearch(radio, display, smiles, group, max, tanimoto) {
  try {
    const queryGroup = group ? `WHERE '${group}' = ANY (array[${taxonomyColumns.join(', ')}])` : '';
    const query = `SELECT DISTINCT structure_smiles FROM data ${queryGroup}`;

    const matchingSmiles = await sendDataPython('public/python/tanimotoSearch.py',query, smiles, tanimoto).catch((error) => {
      console.error("Error while calling Python script:", error);
      return {};
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
  } catch (error) {
    console.error('An error occurred while handling the sim search:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function that finds exact text matches in the database based on the provided criteria.
 * The function searches a specified column for exact matches of the text and returns the data in the desired format.
 * 
 * Parameters:
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - column: The database column where the search should be conducted.
 * - text: The exact text string to be matched against the database column.
 * - max: The maximum number of results to fetch.
 * 
 * Expected outputs:
 * Returns an object containing the data that matches the provided search criteria.
 * 
 * Error Handling:
 * Catches and logs errors related to the search process. Rethrows errors to allow external error handling.
 */async function exactTextMatch(display, column, text, max){
  try {
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
  } catch (error) {
    console.error('An error occurred while executing the exact text match:', error);
    throw error;  // rethrowing the error allows it to be handled further up the call stack if necessary
  }
}


/* async function handleExactTextMatchSPARQL(display, column, text, max){
  switch(display){
    case 'table':
      let query = `
      
      
      `
      return results;
    case 'graph':
      return results;
    default:
      throw new Error(`Unknown display type: ${display}`);
  }
} */

/**
 * Function/Method Comments:
 * An asynchronous function to process and handle exact matches on DBGI data using a SPARQL query.
 * Converts the provided SMILES data to InChI format and formulates a SPARQL query for fetching matching data.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Returns an object containing the results of the SPARQL query.
 * 
 * Error Handling:
 * Catches and logs errors related to the Python script conversion and SPARQL querying. Rethrows errors for further handling.
 */
async function handleSPARQLExactMatch(radio, display, smiles, group, max, tanimoto) {
  try {
    // Convert SMILES to InChI using Python
    let InChi;
    const getInChi = new Promise((resolve, reject) => {
      const process = spawn('python3', ['public/python/smileToInchiKey.py', smiles]);
      
      process.stdout.on('data', (data) => {
        InChi = data.toString();
      });

      process.stderr.on('data', (data) => {
        console.error(`Python stderr: ${data}`);
      });

      process.on('error', (error) => {
        console.error('Error in Python script:', error);
        reject(error);
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
      const graph = await queryGraphDB(query).catch((error) => {
        console.error("Error while querying GraphDB:", error);
        throw error;
      });

      // Render the 'exploreSparql' page with the query results
      return {results: graph.nodes, headers: graph.headers};  
    }
  } catch (error) {
    console.error('An error occurred in handleSPARQLExactMatch:', error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
  return { result: [] };
}

/**
 * Function/Method Comments:
 * An asynchronous function to process and handle substructure searches on DBGI data using a SPARQL query.
 * Queries the database for structures containing the provided substructure.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Returns an object containing the results of the SPARQL query.
 * 
 * Error Handling:
 * Catches and logs errors related to querying GraphDB. Rethrows errors for further handling.
 */
async function handleSPARQLSubSearch(radio, display, smiles, group, max, tanimoto) {
  try {
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
    const graph = await queryGraphDB(query).catch((error) => {
      console.error("Error while querying GraphDB:", error);
      throw error;  // Throw the error to be caught in the outer try-catch block
    });

    console.log(graph);

    // Render the 'exploreSparql' page with the query results
    return {results: graph.nodes, headers: graph.headers};  
  } catch (error) {
    console.error('An error occurred in handleSPARQLSubSearch:', error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
}

/**
 * Function/Method Comments:
 * An asynchronous function that processes and handles similarity searches on DBGI data using a SPARQL query.
 * Constructs the query to find data similar to the provided SMILES string within the specified Tanimoto similarity threshold.
 * 
 * Parameters:
 * - radio: Choice of structure type (either 'structure_inchi' or 'structure_smiles').
 * - display: The format in which data should be displayed ('table' or 'graph').
 * - smiles: The SMILES string to be processed.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * - tanimoto: Threshold for Tanimoto similarity.
 * 
 * Expected outputs:
 * Returns an object containing the results of the SPARQL query.
 * 
 * Error Handling:
 * Catches and logs errors related to the SPARQL query construction and execution. Rethrows errors for external handling.
 */
async function handleSPARQLSimSearch(radio, display, smiles, group, max, tanimoto) {
  try {
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
    const graph = await queryGraphDB(query).catch((error) => {
      console.error("Error while querying GraphDB:", error);
      throw error;  // Throw the error to be caught in the outer try-catch block
    });

    // Render the 'exploreSparql' page with the query results
    return {results: graph.nodes, headers: graph.headers};  
  } catch (error) {
    console.error('An error occurred in handleSPARQLSimSearch:', error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function that fetches data suitable for table display based on the specified criteria.
 * Queries the database using SQL to retrieve rows where the specified column contains the provided structure.
 * 
 * Parameters:
 * - column: The database column where the search should be conducted.
 * - structure: The structure to be matched against the database column.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * 
 * Expected outputs:
 * Returns an object containing the query results.
 * 
 * Error Handling:
 * Catches and logs errors related to SQL query construction and execution. Rethrows errors for external handling.
 */
async function getTableData(column, structure, group, max) {
  try {
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
  } catch (error) {
    console.error(`Error in getTableData function with column ${column}, structure ${structure}, group ${group} and max ${max}:`, error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function that fetches data suitable for graph display based on the specified criteria.
 * Uses SQL to retrieve taxonomical data where the specified column contains the provided structure.
 * 
 * Parameters:
 * - column: The database column for the search.
 * - structure: The structure to be matched against the database column.
 * - group: Classification group to filter by.
 * - max: Maximum number of results to fetch.
 * 
 * Expected outputs:
 * Returns an object containing the graph data and total count of results.
 * 
 * Error Handling:
 * Catches and logs errors related to SQL query construction and execution. Rethrows errors for external handling.
 */
async function getGraphData(column, structure, group, max) {
  try {
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
  } catch (error) {
    console.error(`Error in getGraphData function with column ${column}, structure ${structure}, group ${group} and max ${max}:`, error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
}


/**
 * Function/Method Comments:
 * A function that converts rows of taxonomy data into a hierarchical structure suitable for graph visualization.
 * Iterates over rows and constructs a hierarchical graph structure from taxonomy data.
 * 
 * Parameters:
 * - rows: The rows of taxonomy data to be converted.
 * 
 * Expected outputs:
 * Returns the root of the graph containing the hierarchical structure.
 * 
 * Error Handling:
 * Catches and logs errors encountered during graph construction. Rethrows errors for external handling.
 */
function buildGraphData(rows) {
  try {
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
  } catch (error) {
    console.error("Error in buildGraphData function:", error);
    throw error;  // Re-throw the error if you want to handle it further up the call stack
  }
}


/**
 * Function/Method Comments:
 * An asynchronous function that queries a GraphDB based on the provided SPARQL query and returns the data.
 * Interacts with the GraphDB, passing the query and processing the stream of results to extract the desired data.
 * 
 * Parameters:
 * - query: The SPARQL query to execute against the GraphDB.
 * 
 * Expected outputs:
 * Returns a promise that resolves with the data from the GraphDB.
 * 
 * Error Handling:
 * Catches and logs errors encountered during GraphDB interaction. Rethrows errors for external handling.
 */
async function queryGraphDB(query) {
  try {
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
  } catch (error) {
    console.error('Error while querying GraphDB:', error);
    throw error; // Re-throw the error if you want to handle it further up the call stack
  }
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
          res.render('exploreStructure', { columns, results: results.result, hits: hits  , display: display, source: datasource});
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
