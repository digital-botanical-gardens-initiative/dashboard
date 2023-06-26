process.env.NODE_TLS_REJECT_UNAUTHORIZED = '0'; // should only be used for local development!!! Do not use this method in production environments, as it would expose the app to security risks.

// import modules
const express = require('express');
const dotenv = require("dotenv");
const exploreRoutes = require('./routes/exploreRoutes');
const elementRoutes = require('./routes/elementRoutes');
const organismRoutes = require('./routes/organismRoutes');
const homeRoutes = require('./routes/homeRoutes');
const downloadRoutes = require('./routes/downloadRoutes');



// import dotenv file
dotenv.config();

//express app
const app = express();
const PORT = process.env.PORT || 3000;

// register view engine
app.set('view engine', 'ejs');

//listen for requests
app.listen(PORT, '134.21.20.118', () => console.log('Server listening on 134.21.20.118:' + PORT));

app.use(express.static('public'), exploreRoutes, elementRoutes, organismRoutes, homeRoutes, downloadRoutes);


// Get pages

app.get('/test', (req,res) => {
  res.render('test');
});

app.get('/about', (req,res) => {
    res.render('about', {title: 'About'});
});

app.get('/docu', (req,res) => {
    res.render('docu', {title: 'Documentation'});
});



app.use((req,res) => {
    res.status(404).render('404', {title: '404'});
});

