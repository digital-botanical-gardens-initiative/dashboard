// exploreRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

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
      
      res.render('exploreText', { columns, results: results.rows }); 
    } else {
      res.render('exploreText', { columns });
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
      const smile = req.body.smiles; // Get the SMILES string from the request body
      let results = {}

      if (req.body.exact === 'true') {
        results = await db.query(`SELECT * FROM data WHERE structure_smiles = $1`, [`%${smile}%`]);
      } else {
        results = await db.query(`SELECT * FROM data WHERE structure_smiles ILIKE $1`, [`%${smile}%`]);
      }


      res.render('exploreStructure', { columns, results: results.rows});
    } else {
      res.render('exploreStructure', { columns });
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
