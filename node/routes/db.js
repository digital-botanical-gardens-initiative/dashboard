/**
 * File/Folder Overview:
 * This file sets up a connection pool to a PostgreSQL database using the 'pg' (node-postgres) module.
 * It also provides a utility method for querying the database.
 * 
 * Usage Examples:
 * const db = require('./path_to_this_file');
 * const results = await db.query('SELECT * FROM my_table WHERE id = $1', [1]);
 */

// Import the 'dotenv' module, which loads environment variables from a .env file into process.env
const dotenv = require("dotenv");

// Load the environment variables from the .env file
dotenv.config();

/**
 * Variable Comments:
 * Pool: This is the primary class in the 'pg' module to handle pooling of database connections.
 * Connection pools are used to enhance performance by reusing connection instances.
 */
const Pool = require('pg').Pool;

/**
 * Variable Comments:
 * pool: An instance of Pool that's configured with environment variables.
 * This instance represents our connection pool to the PostgreSQL database.
 * 
 * Configuration details:
 * - host: The PostgreSQL host.
 * - port: The port PostgreSQL runs on.
 * - user: The username for PostgreSQL authentication.
 * - password: The password for PostgreSQL authentication.
 * - database: The name of the PostgreSQL database to connect to.
 * - ssl: A flag to determine if SSL should be used for the connection.
 */
const pool = new Pool({
    host: process.env.PG_HOST,        
    port: process.env.PG_PORT,        
    user: process.env.PG_USER,       
    password: process.env.PG_PASSWORD,
    database: process.env.PG_DATABASE,
    ssl: true,                   
});

/**
 * Function/Method Comments:
 * query: This method provides a utility for querying the PostgreSQL database.
 * Params:
 * - text: The SQL query text.
 * - params: An array of values to be inserted into the SQL query using parameterized query syntax.
 * Expected outputs:
 * - Returns the result of the executed SQL query.
 */
module.exports = {
    query: (text, params) => pool.query(text, params),
};
