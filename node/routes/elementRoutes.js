/**
 * File/Folder Overview:
 * elementRoutes.js -
 * Defines routes and handlers for fetching and displaying elements (chemical structures) and their details.
 */

// Import required modules
const express = require('express');
const router = express.Router();
const db = require('./db');

// Middleware for parsing URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

/**
 * Function/Method Comments:
 * Handles all requests (GET, POST, etc.) to the '/element' route.
 * For POST requests, it processes a search term and fetches matching structure names from the database.
 * For other requests, it just displays structure name examples.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries or other unexpected issues. Sends a 500 status code response for caught errors.
 */
router.all('/element', async (req, res) => {
  try {
    // Example query to retrieve structure names
    const query_example = `
      SELECT structure_nameTraditional
      FROM data
      WHERE structure_nameTraditional IS NOT NULL
      GROUP BY structure_nameTraditional
      ORDER BY COUNT(*) DESC
      LIMIT 12;
    `;

    // Retrieve examples of structure names from the database
    const examples = await db.query(query_example);

    if (req.method === 'POST') {
      // Process POST request
      let searchTerm = req.body.searchTerm.toUpperCase();

      // Query to retrieve distinct structure names matching the search term
      const query = `
        SELECT
        DISTINCT structure_nameTraditional AS molecules
        FROM data
        WHERE structure_nameTraditional ILIKE '${searchTerm}%'
        ORDER BY molecules;
      `;

      // Execute the query and retrieve the result rows
      const { rows } = await db.query(query);

      // Render the 'element' view with the retrieved rows and examples
      res.render('element', { rows, examples: examples.rows });
    } else {
      // Render the 'element' view with only the examples
      res.render('element', { examples: examples.rows });
    }
  } catch (error) {
    console.error('Error fetching elements:', error);
    res.status(500).send('Oops! Looks like something went wrong...');
  }
});

/**
 * Function/Method Comments:
 * Handles GET requests to the '/element/:id' route.
 * Fetches detailed data about an element (chemical structure) based on its ID from the database and displays it.
 * 
 * Parameters:
 * req.params.id - The ID of the element (chemical structure) to fetch details for.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries or other unexpected issues. Sends a 500 status code response for caught errors.
 */
router.get('/element/:id', async (req, res) => {
  try {
    const id = req.params.id;

    // Query to retrieve data for a specific structure name
    const query = `
      SELECT *
      FROM data
      WHERE structure_nameTraditional = '${id}'
      ORDER BY organism_name;
    `;

    // Execute the query and retrieve the result rows
    const { rows } = await db.query(query);

    // Extract relevant data from the first row of the result
    const wikidata = rows[0].structure_wikidata;
    const wikidata_id = wikidata.split('/').pop();
    const formula = rows[0].structure_molecular_formula;
    const smile = rows[0].structure_smiles;
    const smile_2d = rows[0].structure_smiles_2D;
    const kingdom = rows[0].structure_taxonomy_classyfire_01kingdom;
    const superclass = rows[0].structure_taxonomy_classyfire_02superclass;
    const structure_class = rows[0].structure_taxonomy_classyfire_03class;
    const directparent = rows[0].structure_taxonomy_classyfire_04directparent;
    const organisms = [...new Set(rows.map((row) => row.organism_name))];

    // Render the 'elementVisu' view with the retrieved data
    res.render('elementVisu', {
      id,
      wikidata,
      wikidata_id,
      formula,
      smile,
      smile_2d,
      kingdom,
      superclass,
      structure_class,
      directparent,
      organisms,
    });
  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Oops! Looks like something went wrong...');
  }
});

// Export the router module
module.exports = router;
