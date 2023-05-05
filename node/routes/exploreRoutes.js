// exploreRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));




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
    res.send('Error while fetching columns names');
  }
});

router.all('/explore/text', async (req, res) => {
  try {
    const columns = await getTableColumns();

    if (req.method === 'POST') {
      const column = req.body.column;
      const searchTerm = req.body.search;
      const results = await db.query(`SELECT * FROM data WHERE ${column} ILIKE $1`, [`%${searchTerm}%`]);
      res.render('explore_text', { columns, results: results.rows });
    } else {
      res.render('explore_text', { columns });
    }
  } catch (err) {
    console.error(err);
    res.send('Error while fetching columns names');
  }
});


router.get('/explore/structure', (req, res) => {
  res.render('explore_structure');
});



module.exports = router;
