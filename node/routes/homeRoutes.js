// homeRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));

router.get('/', async (req, res) => {
    const result = await db.query(`
        SELECT 
            structure_taxonomy_npclassifier_01pathway as pathway,
            structure_taxonomy_npclassifier_02superclass as superclass,
            structure_taxonomy_npclassifier_03class as class,
            COUNT(DISTINCT structure_nameTraditional) as molecule_count
        FROM data
        GROUP BY pathway, superclass, class
    `);

    const data = result.rows;

    const tree = { name: "Root", children: [] };
    data.forEach(row => {
        let currentLevel = tree.children;
        [row.pathway, row.superclass, row.class].forEach(level => {
            let existingPath = currentLevel.find(d => d.name === level);
            if (existingPath) {
                currentLevel = existingPath.children;
            } else {
                const newPath = { name: level, children: [] };
                currentLevel.push(newPath);
                currentLevel = newPath.children;
            }
        });
    });

    res.render('home', { tree });
});


module.exports = router;

