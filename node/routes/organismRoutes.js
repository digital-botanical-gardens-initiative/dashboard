/**
 * File/Folder Overview:
 * organismRoutes.js - Defines routes and handlers for fetching and displaying organism data.
 * 
 * Usage Examples:
 * To use this module in an Express app:
 * const organismRoutes = require('./path_to_this_file');
 * app.use(organismRoutes);
 */

// Import required modules
const express = require('express');
const router = express.Router();
const db = require('./db');
const axios = require('axios');

// Middleware for parsing URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

/**
 * Function/Method Comments:
 * An asynchronous function to fetch an image of an organism from Wikidata.
 * 
 * Parameters:
 * wikidata_id - The identifier associated with the organism on Wikidata.
 * 
 * Expected outputs:
 * A URL string pointing to the image of the organism or null if an error occurs.
 * 
 * Error Handling:
 * Catches and logs errors related to the Axios request. Returns null if an error is encountered.
 */
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
    return imageUrl;
  } catch (error) {
    console.error(error);
    return null;
  }
}

/**
 * Function/Method Comments:
 * Handles all requests (GET, POST, etc.) to the '/organism' route.
 * For POST requests, it retrieves organisms matching a search term.
 * For other requests, it displays a list of example organisms.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries or other unexpected issues. Sends a 500 status code response for caught errors.
 */
router.all('/organism', async (req, res) => {
  try {
    // Example query to retrieve organism names
    const query_example = `
      SELECT
      organism_name AS organism
      FROM data
      WHERE organism_name IS NOT NULL
      GROUP BY organism
      ORDER BY COUNT(*) DESC
      LIMIT 12;
    `;

    // Retrieve examples of organism names from the database
    const examples = await db.query(query_example);

    if (req.method === 'POST') {
      // Process POST request
      let searchTerm = req.body.searchTerm.toUpperCase();

      // Query to retrieve distinct organism names matching the search term
      const queryOrga = `
        SELECT
        DISTINCT organism_name AS organism
        FROM data
        WHERE organism_name ILIKE '${searchTerm}'
        ORDER BY organism;
      `;

      // Execute the query and retrieve the result rows
      const { rows } = await db.query(queryOrga);

      // Render the 'organism' view with the retrieved rows and examples
      res.render('organism', { rows, examples: examples.rows });
    } else {
      // Render the 'organism' view with only the examples
      res.render('organism', { examples: examples.rows });
    }
  } catch (error) {
    console.error('Error fetching organisms:', error);
    res.status(500).send('Oops! Looks like something went wrong...');
  }
});

/**
 * Function/Method Comments:
 * Handles GET requests to the '/organism/:id' route.
 * Retrieves and displays detailed data for a specific organism based on the provided ID.
 * 
 * Parameters:
 * req.params.id - The ID of the organism to retrieve details for.
 * 
 * Error Handling:
 * Catches and logs errors related to database queries, Axios requests, or other unexpected issues. Sends a 500 status code response for caught errors.
 */
router.get('/organism/:id', async (req, res) => {
  try {
    const id = req.params.id;

    // Query to retrieve data for a specific organism name
    const queryVisu = `
      SELECT *
      FROM data
      WHERE organism_name = '${id}';
    `;

    // Execute the query and retrieve the result rows
    const { rows } = await db.query(queryVisu);

    // Extract relevant data from the first row of the result
    const wikidata = rows[0].organism_wikidata;
    const wikidata_id = wikidata.split('/').pop();
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
    const molecules = [...new Set(rows.map((row) => row.structure_nametraditional))];

    // Fetch the organism image from Wikidata
    const imageUrl = await fetchImage(wikidata_id);

    // Render the 'organismVisu' view with the retrieved data
    res.render('organismVisu', {
      id,
      wikidata,
      wikidata_id,
      domain,
      kingdom,
      phylum,
      organism_class,
      order,
      family,
      tribe,
      genus,
      species,
      varietas,
      molecules,
      imageUrl,
    });
  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Oops! Looks like something went wrong...');
  }
});

// Export the router module
module.exports = router;
