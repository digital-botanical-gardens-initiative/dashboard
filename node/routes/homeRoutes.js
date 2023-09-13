/**
 * homeRoutes.js
 * 
 * This module defines the routes and logic for the home page of the application.
 * It contains functions to interact with a GraphDB, parse the result, and fetch taxon counts.
 * Additionally, it includes caching logic to improve performance on frequent requests.
 */

const express = require('express');
const router = express.Router();
const db = require('./db');

// Importing necessary modules from 'graphdb' package
const {ServerClient, GraphDBServerClient, ServerClientConfig} = require('graphdb').server;
const {RepositoryClientConfig, RDFRepositoryClient} = require('graphdb').repository;
const {RDFMimeType} = require('graphdb').http;
const {SparqlJsonResultParser} = require('graphdb').parser;
const {GetQueryPayload, QueryType} = require('graphdb').query;
const axios = require('axios');

const NodeCache = require( "node-cache" );
const Cache = new NodeCache();

router.use(express.urlencoded({ extended: true }));

/**
 * Queries a GraphDB with a provided SPARQL query and returns the result.
 * 
 * @param {string} query - The SPARQL query to be executed.
 * @returns {Promise} A promise that resolves with the GraphDB data or rejects with an error.
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


/**
 * Parses the raw SPARQL result to extract relevant information.
 * 
 * @param {object} result - Raw result from the GraphDB.
 * @returns {object} An object containing counts of species, genus, families, orders, and kingdoms.
 */

function parseSPARQLResult(result) {
    // Mock-up result parsing
    // Adjust this function to match the actual format of your SPARQL results
    return {
        count_of_species: result[0].count_of_species,
        count_of_genus: result[0].count_of_genus,
        count_of_families: result[0].count_of_families,
        count_of_orders: result[0].count_of_orders,
        count_of_kingdoms: result[0].count_of_kingdoms
    };
}

/**
 * Fetches taxon counts including species, genus, families, orders, and kingdoms.
 * 
 * @returns {Promise} A promise that resolves with the taxon count data or rejects with an error.
 */
async function CountTaxon() {
    try {  
      let query = `PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
      PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
      PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
      PREFIX wdt: <http://www.wikidata.org/prop/direct/>
      PREFIX wd: <http://www.wikidata.org/entity/>
      
      SELECT 
      (COUNT(DISTINCT ?species_name) AS ?count_of_species)
      (COUNT(DISTINCT ?genus_name) AS ?count_of_genus)
      (COUNT(DISTINCT ?family_name) AS ?count_of_families)
      (COUNT(DISTINCT ?order_name) AS ?count_of_orders)
      (COUNT(DISTINCT ?kingdom_name) AS ?count_of_kingdoms)
      WHERE
      {  
          ?material enpkg:has_lab_process ?extract .
          ?material enpkg:has_wd_id ?wd_sp .
          OPTIONAL
          {
                  SERVICE <https://query.wikidata.org/sparql> {
                  ?wd_sp wdt:P225 ?species_name .
                  ?family wdt:P31 wd:Q16521 ;
                      wdt:P105 wd:Q35409 ;
                      wdt:P225 ?family_name ;
                      ^wdt:P171* ?wd_sp .
                  ?genus wdt:P31 wd:Q16521 ;
                      wdt:P105 wd:Q34740 ;
                      wdt:P225 ?genus_name ;
                      ^wdt:P171* ?wd_sp  .
                  ?kingdom wdt:P31 wd:Q16521 ;
                      wdt:P105 wd:Q36732 ;
                      wdt:P225 ?kingdom_name ;
                      ^wdt:P171* ?wd_sp .
                  ?order wdt:P31 wd:Q16521 ;
                      wdt:P105 wd:Q36602 ;
                      wdt:P225 ?order_name ;
                      ^wdt:P171* ?wd_sp .
              }
          }
      } 
      `;
  
        // Query the GraphDB using the SPARQL query and get the results
        const result = await queryGraphDB(query).catch((error) => {
            console.error("Error while querying GraphDB:", error);
            throw error;  // Throw the error to be caught in the outer try-catch block
        });

        // Parse the SPARQL result and return it
        return parseSPARQLResult(result.nodes);
    } catch (error) {
      console.error('An error occurred:', error);
      throw error;  // Re-throw the error if you want to handle it further up the call stack
    }
  }

  

/**
 * Fetches the species count from a specified project using the iNaturalist API.
 * 
 * @async
 * @param {string} projectSlug - The unique identifier (slug) for the project in the iNaturalist platform.
 * @returns {number|null} Returns the total species count for the given project. If an error occurs, it returns null.
 * @throws Will log an error to the console if the API call or data extraction fails.
 */
async function getSpeciesCountFromProject(projectSlug) {
  const baseURL = 'https://api.inaturalist.org/v1/observations/species_counts';
  
  try {
      const response = await axios.get(baseURL, {
          params: {
              project_id: projectSlug,
              verifiable: false
          }
      });
      
      if (response.data && response.data.results) {
          return response.data.total_results;
      } else {
          throw new Error('Failed to fetch data.');
      }
  } catch (error) {
      console.error('Error fetching species count:', error);
      return null;
  }
}

const projectSlug = 'digital-botanical-gardens-initiative';
getSpeciesCountFromProject(projectSlug).then(count => {
  console.log(`Number of species in project ${projectSlug}:`, count);
});


/**
 * Router handler for the home route.
 * 
 * Checks the cache for data first, if not found, then queries the database.
 * Sends the resulting data to the 'home' view for rendering.
 */
router.all('/', async (req, res) => {
    try {
        // Fetch the data from CountTaxon
        let data = Cache.get( "myKey" );  // Initialize data with cached value
        if (data == undefined){    
          data = await CountTaxon();  // Assign directly to the data variable
          Cache.set("myKey", data, 100000);
        }
        // Assume data returns a row of results as an object e.g., { species_count: 10, order_count: 20, ...}
        // Pass this data into the 'home' view
        res.render('home', {results: data});
    } catch (error) {
        console.error('An error occurred:', error);
        res.status(500).send('An error occurred');
    }
});



module.exports = router;

