// exploreRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

const getTableColumns = async () => {
  const tableColumns = await db.query(`
    SELECT column_name
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'data'
  `);
  return tableColumns.rows;
};

router.get('/explore', async (req, res) => {
  try {
    const columns = await getTableColumns();
    res.render('explore', { columns });
  } catch (err) {
    console.error(err);
    res.send('Error');
  }
});

router.post('/explore', async (req, res) => {
  try {
    //const column = req.body.column;
    console.log(req.body);
    const searchTerm = req.body.search;
    const results = await db.query(`SELECT * FROM data WHERE ${column} ILIKE $1`, [`%${searchTerm}%`]);
    res.render('explore', { results: results.rows });
  } catch (err) {
    console.error(err);
    res.send('Error');
  }
});


module.exports = router;
