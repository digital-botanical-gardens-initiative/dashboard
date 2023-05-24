// elementRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));

router.all('/element', async (req, res) => {
  try {

    const query_example = `
    SELECT structure_nameTraditional
    FROM data
    WHERE structure_nameTraditional IS NOT NULL
    GROUP BY structure_nameTraditional
    ORDER BY COUNT(*) DESC
    LIMIT 12;
    `;

    const examples  = await db.query(query_example);

    if (req.method === 'POST') {
      let searchTerm = req.body.searchTerm.toUpperCase();
      const query = `
          SELECT
          DISTINCT structure_nameTraditional AS molecules
          FROM data
          WHERE structure_nameTraditional ILIKE '${searchTerm}%'
          ORDER BY molecules;
      `;

      const { rows } = await db.query(query);

      res.render('element', { rows , examples: examples.rows});
    } else {
      res.render('element', {examples : examples.rows});
    }
  } catch (error) {
      console.error('Error fetching elements:', error);
      res.status(500).send('Oops! Looks like something went wrong...')
  }
});
  
  
router.get('/element/:id', async (req,res) =>{
    try{
      const id = req.params.id;
      const query = `
        SELECT *
        FROM data
        WHERE structure_nameTraditional = '${id}'
        ORDER BY organism_name;`;
  
      const { rows } = await db.query(query);
      const wikidata = rows[0].structure_wikidata;
      const wikidata_id = wikidata.split("/").pop();
      const formula = rows[0].structure_molecular_formula;
      const smile = rows[0].structure_smiles;
      const smile_2d = rows[0].structure_smiles_2D;
      const kingdom = rows[0].structure_taxonomy_classyfire_01kingdom;
      const superclass = rows[0].structure_taxonomy_classyfire_02superclass;
      const structure_class = rows[0].structure_taxonomy_classyfire_03class;
      const directparent = rows[0].structure_taxonomy_classyfire_04directparent;
      const organisms = [...new Set(rows.map(row => row.organism_name))];    
  
    
    res.render('elementVisu', { id, wikidata, wikidata_id, formula, smile, smile_2d, kingdom, superclass, structure_class, directparent, organisms });
  
    } catch (error) {
      console.error('Error:', error);
      res.status(500).send('Oops! Looks like something went wrong...')
    }
  })

  module.exports = router;
