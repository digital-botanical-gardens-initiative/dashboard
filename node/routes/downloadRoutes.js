// downloadRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');
const fs = require('fs');
const json2csv = require('json2csv').parse;
const jsonfile = require('jsonfile');
/* const pgDump = require('pg-dump');
 */
router.use(express.urlencoded({ extended: true }));

router.get('/download', (req, res) => {
  try {
      res.render('download', {title: 'Download'});
  } catch (error) {
      console.error(error);
      res.status(500).send("Server Error");
  }
});

router.get('/download/csv', async (req, res) => {
  try {
      const data = await db.query('SELECT * FROM data');
      const csv = json2csv(data.rows);
      fs.writeFileSync('./dbgi.csv', csv);
      res.download('./dbgi.csv');
  } catch (error) {
      console.error(error);
      if(error.code === 'ENOENT') {
          res.status(404).send("File not found");
      } else {
          res.status(500);
          res.redirect('/download');
      }
  }
});

router.get('/download/json', async (req, res) => {
  try {
      const data = await db.query('SELECT * FROM data');
      jsonfile.writeFileSync('./dbgi.json', data.rows);
      res.download('./dbgi.json');
  } catch (error) {
      console.error(error);
      if(error.code === 'ENOENT') {
          res.status(404).send("File not found");
      } else {
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

module.exports = router;


