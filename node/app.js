/**
 * File/Folder Overview:
 * This file sets up and configures the Express application server for our project.
 * It handles initialization of middlewares, route handling, error pages, and server startup.
 */

// import modules
const express = require('express');

const dotenv = require("dotenv");
const exploreRoutes = require('./routes/exploreRoutes');
const elementRoutes = require('./routes/elementRoutes');
const organismRoutes = require('./routes/organismRoutes');
const homeRoutes = require('./routes/homeRoutes');
const downloadRoutes = require('./routes/downloadRoutes');
const KGRoutes = require('./routes/KGRoutes');



// Load environment variables from the .env file
dotenv.config();

// Initialize the Express application
const app = express();
const PORT = process.env.PORT || 3000;

// Setup the view engine for the application
app.set('view engine', 'ejs');

// Start the server and listen for incoming requests
// Function/Method Comments: This function starts the server on the specified IP and PORT.
app.listen(PORT, '134.21.20.118', () => console.log('Server listening on 134.21.20.118:' + PORT));

// Middleware setup for serving static files and routing
app.use(express.static('public'), exploreRoutes, elementRoutes, organismRoutes, homeRoutes, downloadRoutes, KGRoutes);

/**
 * Function/Method Comments:
 * The function below renders an 'About' page when the '/about' route is hit.
 * Params:
 * - req: The request object containing request information.
 * - res: The response object used to send back the required response.
 * Expected outputs:
 * - Renders the 'about' view.
 */
app.get('/about', (req,res) => {
    res.render('about', {title: 'About'});
});

/**
 * Function/Method Comments:
 * The function below renders a 'Documentation' page when the '/docu' route is hit.
 * Params:
 * - req: The request object containing request information.
 * - res: The response object used to send back the required response.
 * Expected outputs:
 * - Renders the 'docu' view.
 */
app.get('/docu', (req,res) => {
    res.render('docu', {title: 'Documentation'});
});

/**
 * Error Handling:
 * If none of the routes match, this function sends a 404 page.
 */
app.use((req,res) => {
    res.status(404).render('404', {title: '404'});
});

