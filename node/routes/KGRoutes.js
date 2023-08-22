/**
 * File/Folder Overview:
 * KGRoutes.js - Defines routes and handlers for querying a GraphDB and rendering knowledge graph (KG) views.
 * 
 * Usage Examples:
 * To use this module in an Express app:
 * const kgRoutes = require('./path_to_this_file');
 * app.use(kgRoutes);
 */

// Import required modules
const express = require('express');
const router = express.Router();
const dotenv = require("dotenv");
dotenv.config();
const { ServerClient, GraphDBServerClient, ServerClientConfig } = require('graphdb').server;
const { RepositoryClientConfig, RDFRepositoryClient } = require('graphdb').repository;
const { RDFMimeType } = require('graphdb').http;
const { SparqlJsonResultParser } = require('graphdb').parser;
const { GetQueryPayload, QueryType } = require('graphdb').query;

// Middleware for parsing URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

/**
 * Function/Method Comments:
 * An asynchronous function to send queries to GraphDB and retrieve results.
 * 
 * Parameters:
 * query - A SPARQL query string.
 * 
 * Expected outputs:
 * An object with nodes (graph data) and headers (bindings).
 * 
 * Side effects:
 * Sends a query to the GraphDB and retrieves results. 
 * 
 * Error Handling:
 * If there's an error in the stream or query process, the promise is rejected with the respective error.
 */
async function queryGraphDB(query) {
  const endpoint = process.env.ENDPOINT_GRAPHDB;
  const readTimeout = 30000;
  const writeTimeout = 30000;
  const config = new RepositoryClientConfig(endpoint)
    .setEndpoints([process.env.REPO_GRAPHDB])
    .setHeaders({
      'Accept': RDFMimeType.TURTLE
    })
    .setReadTimeout(readTimeout)
    .setWriteTimeout(writeTimeout);

  const repository = new RDFRepositoryClient(config);

  repository.registerParser(new SparqlJsonResultParser());

  const payload = new GetQueryPayload()
    .setQuery(query)
    .setQueryType(QueryType.SELECT)
    .setResponseType(RDFMimeType.SPARQL_RESULTS_JSON)
    .setLimit(100);

  return new Promise((resolve, reject) => {
    const nodes = [];
    let headers = [];

    repository.query(payload).then((stream) => {
      stream.on('data', (bindings) => {
        const node = {};
        for (const binding in bindings) {
          if (bindings.hasOwnProperty(binding)) {
            node[binding] = bindings[binding].value; // use binding as id and bindings[binding].value as label
            if (!headers.includes(binding)) {
              headers.push(binding); // add the binding to the headers if it's not already there
            }
          }
        }
        nodes.push(node);
      });
      stream.on('end', () => {
        // Resolve the promise with the graph data and the headers when the stream ends
        resolve({ nodes: nodes, headers: headers });
      });

      stream.on('error', (error) => {
        // Reject the promise if there's an error
        reject(error);
      });
    });
  });
}

/**
 * Function/Method Comments:
 * Handles all requests (GET, POST, etc.) to the '/KG' route.
 * For POST requests, retrieves a SPARQL query from the request body and fetches corresponding data from GraphDB.
 * For other requests, displays the default 'KG' view.
 * 
 * Error Handling:
 * Catches and logs errors related to GraphDB queries or other unexpected issues. Sends a 500 status code response for caught errors.
 */
router.all('/KG', async (req, res) => {
  try {
    if (req.method === 'POST') {
      // Get the query from the request body
      const query = req.body.query;

      // Query the GraphDB and retrieve the graph data
      const graph = await queryGraphDB(query);

      // Render the 'KG' view with the retrieved graph data
      res.render('KG', { results: graph.nodes, headers: graph.headers });
    } else {
      // Render the 'KG' view with the title
      res.render('KG', { title: 'KG' });
    }
  } catch (error) {
    console.error(error);
    res.status(500).send("Server Error");
  }
});

// Export the router module
module.exports = router;
