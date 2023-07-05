// downloadRoutes.js

// Import required modules
const express = require('express');
const router = express.Router();
const db = require('./db');
const fs = require('fs');
const json2csv = require('json2csv').parse;
const jsonfile = require('jsonfile');

// Middleware for parsing URL-encoded bodies
router.use(express.urlencoded({ extended: true }));

// Handler for the '/download' route
router.get('/download', (req, res) => {
  try {
    // Render the 'download' view with the provided title
    res.render('download', { title: 'Download' });
  } catch (error) {
    console.error(error);
    res.status(500).send("Server Error");
  }
});

// Handler for the '/download/csv' route
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

// Handler for the '/download/json' route
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


