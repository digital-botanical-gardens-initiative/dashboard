const express = require('express');
const dotenv = require("dotenv");
const RDKit=require('rdkit'); 

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

app.use(express.static('public'));

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
      ARRAY_AGG(structure_nameTraditional) AS molecules
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
    res.status(500).send('Error fetching elements');
  }
});


app.get('/element/:id', async (req,res) =>{
  try{
    const id = req.params.id;
    const query = `
      SELECT *
      FROM data
      WHERE structure_nameTraditional = '${id}';`;

    const { rows } = await pool.query(query);
    const wikidata = rows[0].structure_wikidata;
    const wikidata_id = wikidata.split("/").pop();
    const formula = rows[0].structure_molecular_formula;
    const smile = rows[0].structure_smiles;
    const smile_2d = rows[0].structure_smiles_2D;
  
  res.render('element_visu', { id, wikidata, wikidata_id, formula, smile, smile_2d });

  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Error in the element visualization')
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
      ARRAY_AGG(organism_name) AS organisms
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
    res.status(500).send('Error fetching organisms');
  }
});

app.get('/organism/<id>', async (req,res) =>{
  try{

  } catch (error) {
    console.error('Error:', error);
    res.status(500).send('Error in the organism visualization')
  }
})

app.get('/explore', (req,res) => {
    res.render('explore', {title: 'Explore'});
});

app.use((req,res) => {
    res.status(404).render('404', {title: '404'});
});
