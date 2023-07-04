// downloadRoutes.js
const express = require('express');
const router = express.Router();
const dotenv = require("dotenv");
const { NONE } = require('graphdb/lib/transaction/transaction-isolation-level');
dotenv.config();
const {ServerClient, GraphDBServerClient, ServerClientConfig} = require('graphdb').server;
const {RepositoryClientConfig, RDFRepositoryClient} = require('graphdb').repository;
const {RDFMimeType} = require('graphdb').http;
const {SparqlJsonResultParser} = require('graphdb').parser;
const {GetQueryPayload, QueryType} = require('graphdb').query;


router.use(express.urlencoded({ extended: true }));

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
    const edges = [];

    repository.query(payload).then((stream) => {
      stream.on('data', (bindings) => {
        console.log(bindings);
        for (const binding in bindings) {
          console.log(binding);
          if (bindings.hasOwnProperty(binding)) {
            nodes.push({id: bindings[binding].value, label: bindings[binding].value});
          }
        }
      });
      

      stream.on('end', () => {
        // Resolve the promise with the graph data when the stream ends
        resolve({nodes: nodes, edges: edges});
      });

      stream.on('error', (error) => {
        // Reject the promise if there's an error
        reject(error);
      });
    });
  });
}

// downloadRoutes.js
router.all('/KG', async (req, res) => {
  try {
    if (req.method === 'POST') {
      const query = req.body.query;
      const graph = await queryGraphDB(query);

      res.render('KG', {results: graph.nodes});  // pass nodes array to the view
    } else {
      res.render('KG', {title: 'KG'});
    }
  } catch (error) {
    console.error(error);
    res.status(500).send("Server Error");
  }
});



module.exports = router;

