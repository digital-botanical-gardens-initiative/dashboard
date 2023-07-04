// downloadRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');
const {ServerClient, GraphDBServerClient, ServerClientConfig} = require('graphdb').server;
const {RepositoryClientConfig, RDFRepositoryClient} = require('graphdb').repository;
const {RDFMimeType} = require('graphdb').http;
const {SparqlJsonResultParser} = require('graphdb').parser;
const {GetQueryPayload, QueryType} = require('graphdb').query;

const dotenv = require("dotenv");
dotenv.config();
require('dotenv').config();


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
                        .setResponseType(RDFMimeType.SPARQL_RESULTS_XML)
                        .setLimit(100);

  return repository.query(payload).then((stream) => {
    stream.on('data', (bindings) => {
      console.log(bindings);
    });
      stream.on('end', () => {
      // handle end of the stream
    });
  });
}


router.all('/KG', (req, res) => {
  try {
    if (req.method === 'POST') {
      const query = req.body.query;
      queryGraphDB(query);
    }

    res.render('KG', {title: 'KG'});
  } catch (error) {
      console.error(error);
      res.redirect('/KG');
      res.status(500).send("Server Error");
  }

}); 

module.exports = router;

