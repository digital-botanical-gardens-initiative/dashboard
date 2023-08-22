/**
 * File/Folder Overview:
 * downloadRoutes.js
 * This file contains routes and handlers for downloading database content in various formats, such as CSV and JSON.
 */

// Import required modules
const express = require('express');
const router = express.Router();
const db = require('./db');
const fs = require('fs');
const json2csv = require('json2csv').parse;
const jsonfile = require('jsonfile');

// Middleware for parsing URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

/**
 * Function/Method Comments:
 * Handles GET requests to the '/download' route.
 * Renders the 'download' view .
 * Error Handling:
 * Catches and logs errors, and sends a 500 status code response in case of unexpected errors.
 */
router.get('/download', (req, res) => {
  try {
    // Render the 'download' view with the provided title
    res.render('download', { title: 'Download' });
  } catch (error) {
    console.error(error);
    res.status(500).send("Server Error");
  }
});

/**
 * Function/Method Comments:
 * Handles GET requests to the '/download/csv' route.
 * Fetches data from the database, converts it to CSV format, writes it to a file, and initiates a download for the user.
 * Error Handling:
 * Catches and logs errors. Specific handling for file not found (ENOENT) errors, and other errors redirect to '/download' with a 500 status code.
 */
router.get('/download/csv', async (req, res) => {
  try {
    // Query the database to retrieve data
    const data = await db.query('SELECT * FROM data');

    // Convert the data to CSV format
    const csv = json2csv(data.rows);

    // Write the CSV data to a file named 'dbgi.csv'
    fs.writeFileSync('./dbgi.csv', csv);

    // Download the 'dbgi.csv' file
    res.download('./dbgi.csv');
  } catch (error) {
    console.error(error);
    if (error.code === 'ENOENT') {
      // If the file is not found, send a 404 error response
      res.status(404).send("File not found");
    } else {
      // If there is a server error, send a 500 error response and redirect to '/download'
      res.status(500);
      res.redirect('/download');
    }
  }
});

/**
 * Function/Method Comments:
 * Handles GET requests to the '/download/json' route.
 * Fetches data from the database, writes it to a JSON file, and initiates a download for the user.
 * Error Handling:
 * Catches and logs errors. Specific handling for file not found (ENOENT) errors, and other errors redirect to '/download' with a 500 status code.
 */
router.get('/download/json', async (req, res) => {
  try {
    // Query the database to retrieve data
    const data = await db.query('SELECT * FROM data');

    // Write the JSON data to a file named 'dbgi.json'
    jsonfile.writeFileSync('./dbgi.json', data.rows);

    // Download the 'dbgi.json' file
    res.download('./dbgi.json');
  } catch (error) {
    console.error(error);
    if (error.code === 'ENOENT') {
      // If the file is not found, send a 404 error response
      res.status(404).send("File not found");
    } else {
      // If there is a server error, send a 500 error response and redirect to '/download'
      res.status(500);
      res.redirect('/download');
    }
  }
});


/* router.get('/download/sql', async (req, res) => {
  const dump = await pgDump.dump({
    user: process.env.PG_USER,
    password: process.env.PG_PASSWORD,
    host: process.env.PG_HOST,
    port: process.env.PG_PORT,
    database: process.env.PG_DATABASE
  });

  fs.writeFileSync('./dbgi.sql', dump);
  res.download('./dbgi.sql');
}); */

// Export the router module
module.exports = router;


