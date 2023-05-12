process.env.NODE_TLS_REJECT_UNAUTHORIZED = '0'; // should only be used for local development!!! Do not use this method in production environments, as it would expose the app to security risks.

const express = require('express');
const dotenv = require("dotenv");
const RDKit=require('rdkit'); 
const exploreRoutes = require('./routes/exploreRoutes');
const elementRoutes = require('./routes/elementRoutes');
const organismRoutes = require('./routes/organismRoutes');


// import dotenv file
dotenv.config();

//express app
const app = express();

// register view engine
app.set('view engine', 'ejs');

//listen for requests
app.listen(3000, '134.21.20.118', () => console.log('Server listening on 134.21.20.118:3000'));

app.use(express.static('public'), exploreRoutes, elementRoutes, organismRoutes);

/* app.use(express.urlencoded({ extended: true }));
app.use(express.json());


require('dotenv').config();

const Pool = require('pg').Pool
const pool = new Pool({
    host: process.env.PG_HOST,
    port: process.env.PG_PORT,
    user: process.env.PG_USER,
    password: process.env.PG_PASSWORD,
    database: process.env.PG_DATABASE,
    ssl: true,
}) */

/* function search(req, res, next) {
  // user's search term
  var searchTerm = req.query.search;
  var column = req.query.column;

    let query = 'SELECT * FROM data';

    if (searchTerm != '' && column != ''){
      query = `SELECT * FROM data WHERE ` + column + `= '` + searchTerm + `'`; 
    }

  pool.query(query, (err, result) => {
    if(err){
      req.searchResult = '';
      req.searchTerm = '';
      req.column = '';
      next();
    }

    req.searchResult = result;
    req.searchTerm = searchTerm;
    req.column = column;
    
    next();
  })
} */


app.get('/', (req,res) => {
    res.render('home', { title: 'Home'});
});

app.get('/test', (req,res) => {
  res.render('test');
});

app.get('/about', (req,res) => {
    res.render('about', {title: 'About'});
});

app.get('/docu', (req,res) => {
    res.render('docu', {title: 'Documentation'});
});

/* app.get('/api/molecules', async (req, res) => {
  try {
    const query = `SELECT structure_taxonomy_classyfire_01kingdom, 
                          structure_taxonomy_classyfire_02superclass, 
                          structure_taxonomy_classyfire_03class, 
                          structure_taxonomy_classyfire_04directparent, 
                          structure_nameTraditional 
                          FROM data`;
    const { rows } = await pool.query(query);
    res.json(rows);
  } catch (error) {
    console.error('Error fetching data:', error);
    res.status(500).json({ error: 'Internal server error' });
  }
}); */


app.use((req,res) => {
    res.status(404).render('404', {title: '404'});
});
