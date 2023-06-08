// organismRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');
const axios = require('axios');


router.use(express.urlencoded({ extended: true }));


async function fetchImage(wikidata_id) {
  const query_wikidata = encodeURIComponent(`
    SELECT ?image WHERE {
      wd:${wikidata_id} wdt:P18 ?image.
    }
  `);

  const url = `https://query.wikidata.org/sparql?format=json&query=${query_wikidata}`;

  try {
    const response = await axios.get(url);
    const imageUrl = response.data.results.bindings[0].image.value;
  } catch (error) {
    console.error(error);
  }
}


router.all('/organism', async (req, res) => {
    try {
      const query_example = `
        SELECT
        organism_name AS organism
        FROM data
        WHERE organism_name IS NOT NULL
        GROUP BY organism
        ORDER BY COUNT(*) DESC
        LIMIT 12;
      `;

      const examples = await db.query(query_example);
  
      if (req.method == 'POST'){
        let searchTerm = req.body.searchTerm.toUpperCase();

        const queryOrga = `
          SELECT
          DISTINCT organism_name AS organism
          FROM data
          WHERE organism_name ILIKE '${searchTerm}'
          ORDER BY organism;
        `
        const { rows } = await db.query(queryOrga);
        
  
        res.render('organism', { rows, examples: examples.rows });
      } else {
        res.render('organism', {examples: examples.rows});
      }
    } catch (error) {
      console.error('Error fetching organisms:', error);
      res.status(500).send('Oops! Looks like something went wrong...')
    }
  });

router.get('/organism/:id', async (req,res) =>{
    try{
      const id = req.params.id;
      const queryVisu = `
        SELECT *
        FROM data
        WHERE organism_name = '${id}';`;
  
      const { rows } = await db.query(queryVisu);
      const wikidata = rows[0].organism_wikidata;
      const wikidata_id = wikidata.split("/").pop();
      const domain = rows[0].organism_taxonomy_01domain;
      const kingdom = rows[0].organism_taxonomy_02kingdom;
      const phylum = rows[0].organism_taxonomy_03phylum;
      const organism_class = rows[0].organism_taxonomy_04class;
      const order = rows[0].organism_taxonomy_05order;
      const family = rows[0].organism_taxonomy_06family;
      const tribe = rows[0].organism_taxonomy_07tribe;
      const genus = rows[0].organism_taxonomy_08genus;
      const species = rows[0].organism_taxonomy_09species;
      const varietas = rows[0].organism_taxonomy_10varietas;
      const molecules = [...new Set(rows.map(row => row.structure_nametraditional))];  


      const query_wikidata = encodeURIComponent(`
      SELECT ?image WHERE {
        wd:${wikidata_id} wdt:P18 ?image.
      }
    `);
  
      const url = `https://query.wikidata.org/sparql?format=json&query=${query_wikidata}`;

      const response = await axios.get(url);
      var imageUrl = null;

      if (response.data.results.bindings.length !== 0) {
        imageUrl = response.data.results.bindings[0].image.value;
      }

    
    res.render('organismVisu', {id, wikidata, wikidata_id, 
                                domain, kingdom, phylum, organism_class, 
                                order, family, tribe, genus, 
                                species, varietas, molecules, imageUrl});
  
    } catch (error) {
      console.error('Error:', error);
      res.status(500).send('Oops! Looks like something went wrong...')
      }
  });


  module.exports = router;