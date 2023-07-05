// Import the 'dotenv' module, which loads environment variables from a .env file into process.env
const dotenv = require("dotenv");

// Load the environment variables from the .env file
dotenv.config();

// Create a new connection pool using the 'pg' (node-postgres) module
const Pool = require('pg').Pool;

// Configure the connection pool using environment variables
const pool = new Pool({
    host: process.env.PG_HOST,        // PostgreSQL host
    port: process.env.PG_PORT,        // PostgreSQL port
    user: process.env.PG_USER,        // PostgreSQL username
    password: process.env.PG_PASSWORD,// PostgreSQL password
    database: process.env.PG_DATABASE,// PostgreSQL database name
    ssl: true,                        // Use SSL for the connection
});

// Export an object with a 'query' method, which allows querying the PostgreSQL database
module.exports = {
    query: (text, params) => pool.query(text, params),
};
