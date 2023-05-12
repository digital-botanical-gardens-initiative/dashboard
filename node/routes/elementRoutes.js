// elementRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));

router.get('/element', async (req, res) => {
    try {
      const query = `
        SELECT
        UPPER(SUBSTRING(structure_nameTraditional, 1, 1)) AS first_letter,
        ARRAY_AGG(DISTINCT structure_nameTraditional) AS molecules
        FROM data
        WHERE structure_nameTraditional IS NOT NULL
        GROUP BY first_letter
        ORDER BY first_letter;
      `;
  
      const { rows } = await db.query(query);
  
      // Convert the result to an object
      const groupedMolecules = {};
      rows.forEach(row => {
        groupedMolecules[row.first_letter] = row.molecules;
      });
  
      res.render('element', { groupedMolecules });
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
