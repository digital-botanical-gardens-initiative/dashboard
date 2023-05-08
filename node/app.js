const express = require('express');
const dotenv = require("dotenv");
const RDKit=require('rdkit'); 
const exploreRoutes = require('./routes/exploreRoutes');

// Import the initRDKit function from the external script file
const { initRDKit } = require("./routes/moleculeVisuRoutes.js");


// import dotenv file
dotenv.config();

//express app
const app = express();

// register view engine
app.set('view engine', 'ejs');

//listen for requests
app.listen(3000, '134.21.20.118', () => console.log('Server listening on 134.21.20.118:3000'));

app.use(express.static('public'), exploreRoutes);

app.use(express.urlencoded({ extended: true }));
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
})

function search(req, res, next) {
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
}


app.get('/', (req,res) => {
    res.render('home', { title: 'Home'});
});

app.get('/about', (req,res) => {
    res.render('about', {title: 'About'});
});

app.get('/docu', (req,res) => {
    res.render('docu', {title: 'Documentation'});
});

app.get('/element', async (req, res) => {
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

    const { rows } = await pool.query(query);

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


app.get('/element/:id', async (req,res) =>{
  try{
    const id = req.params.id;
    const query = `
      SELECT *
      FROM data
      WHERE structure_nameTraditional = '${id}'
      ORDER BY organism_name;`;

    const { rows } = await pool.query(query);
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

  
  res.render('element_visu', { id, wikidata, wikidata_id, formula, smile, smile_2d, kingdom, superclass, structure_class, directparent, organisms });

  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Oops! Looks like something went wrong...')
  }
})

app.post('/element/:id', function(req, res){
  var mol = RDKit.Molecule.fromSmiles( smiles );
  var molwt = mol.getMW();
  var mol2d = mol.Drawing2D();
  var remol = mol2d.replace( /svg:/g, ''  );

  res.send( "molwt is:"+molwt+"<br><br>"+"smiles is:"+ smiles+"<br><br>"+remol );
});



app.get('/organism', async (req, res) => {
  try {
    const query = `
      SELECT
      UPPER(SUBSTRING(organism_name, 1, 1)) AS first_letter,
      ARRAY_AGG(DISTINCT organism_name) AS organisms
      FROM data
      WHERE organism_name IS NOT NULL
      GROUP BY first_letter
      ORDER BY first_letter;
    `;

    const { rows } = await pool.query(query);

    // Convert the result to an object
    const groupedOrganisms = {};
    rows.forEach(row => {
      groupedOrganisms[row.first_letter] = row.organisms;
    });

    res.render('organism', { groupedOrganisms });
  } catch (error) {
    console.error('Error fetching organisms:', error);
    res.status(500).send('Oops! Looks like something went wrong...')
  }
});

app.get('/organism/:id', async (req,res) =>{
  try{
    const id = req.params.id;
    const query = `
      SELECT *
      FROM data
      WHERE organism_name = '${id}';`;

    const { rows } = await pool.query(query);
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

  
  res.render('organism_visu', { id, wikidata, wikidata_id, 
                                domain, kingdom, phylum, organism_class, 
                                order, family, tribe, genus, 
                                species, varietas, molecules});

  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Oops! Looks like something went wrong...')
    }
})


app.get('/explore', (req,res) => {
  try{
    var searchResult = req.searchResult;

    res.render('explore', {title: 'Explore',
                          results: searchResult.length,
                          searchTerm: req.searchTerm,
                          searchResult: searchResult,
                          column: req.column});

  } catch (error) {
    console.error('Error: ', error);
    res.status(500).send('Oops! Looks like something went wrong...')
  }
});

app.get('/explore/text', (req,res) => {
  try{
    res.render('explore_text');

  } catch (error) {
    console.error('Error: ', error);
    res.status(500).send('Oops! Looks like something went wrong...')
  }
});

app.get('/api/molecules', async (req, res) => {
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
});



app.use((req,res) => {
    res.status(404).render('404', {title: '404'});
});
